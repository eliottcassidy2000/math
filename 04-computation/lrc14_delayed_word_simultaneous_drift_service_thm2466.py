#!/usr/bin/env python3
"""Independent exact referee for THM-2466.

The proof is analytic; this dependency-free referee checks the two-BV bound
on a rational step-function bank, the simultaneous clock threshold,
the THM-2459 invoice, and the small-clock/zero-mean hostile boundaries.
"""

from fractions import Fraction as F


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def mean_step(segments: tuple[tuple[F, F, F], ...]) -> F:
    return sum((right - left) * value for left, right, value in segments)


def validate_step(segments: tuple[tuple[F, F, F], ...]) -> None:
    require(bool(segments), "step function needs a segment")
    require(segments[0][0] == 0 and segments[-1][1] == 1, "step endpoints")
    for index, (left, right, _) in enumerate(segments):
        require(left < right, "positive step interval")
        if index:
            require(segments[index - 1][1] == left, "step partition gap")


def variation_step(segments: tuple[tuple[F, F, F], ...]) -> F:
    validate_step(segments)
    values = [value for _, _, value in segments]
    return sum(abs(values[index] - values[index - 1]) for index in range(len(values)))


def validate_word(intervals: tuple[tuple[F, F], ...]) -> None:
    previous = F(0)
    for left, right in intervals:
        require(F(0) <= left < right <= F(1), "word interval")
        require(left >= previous, "overlapping word intervals")
        previous = right


def word_mean(intervals: tuple[tuple[F, F], ...]) -> F:
    validate_word(intervals)
    return sum(right - left for left, right in intervals)


def word_value(intervals: tuple[tuple[F, F], ...], point: F) -> int:
    require(F(0) <= point < F(1), "word evaluation point")
    return int(any(left <= point < right for left, right in intervals))


def word_variation(intervals: tuple[tuple[F, F], ...]) -> F:
    validate_word(intervals)
    breaks = sorted({F(0), F(1), *(x for interval in intervals for x in interval)})
    values = [
        word_value(intervals, (left + right) / 2)
        for left, right in zip(breaks, breaks[1:])
    ]
    return F(sum(abs(values[index] - values[index - 1]) for index in range(len(values))))


def integrate_pulled_word(
    segments: tuple[tuple[F, F, F], ...],
    intervals: tuple[tuple[F, F], ...],
    frequency: int,
) -> F:
    """Integral f(x) 1_W(frequency*x) dx by exact interval intersection."""

    validate_step(segments)
    validate_word(intervals)
    require(frequency >= 1, "positive frequency")
    total = F(0)
    for left, right, value in segments:
        for cell in range(frequency):
            for word_left, word_right in intervals:
                pull_left = (cell + word_left) / frequency
                pull_right = (cell + word_right) / frequency
                overlap = min(right, pull_right) - max(left, pull_left)
                if overlap > 0:
                    total += overlap * value
    return total


def covariance_bound(
    segments: tuple[tuple[F, F, F], ...],
    intervals: tuple[tuple[F, F], ...],
    frequency: int,
) -> tuple[F, F]:
    alpha = word_mean(intervals)
    actual = integrate_pulled_word(segments, intervals, frequency)
    error = abs(actual - alpha * mean_step(segments))
    bound = word_variation(intervals) * variation_step(segments) / (12 * frequency)
    return error, bound


def main() -> None:
    words = (
        ((F(0), F(1, 3)),),
        ((F(1, 7), F(2, 7)), (F(4, 7), F(6, 7))),
        ((F(0), F(1)),),
    )
    functions = (
        ((F(0), F(1), F(1)),),
        ((F(0), F(1, 2), F(3)), (F(1, 2), F(1), F(-1))),
        (
            (F(0), F(1, 5), F(0)),
            (F(1, 5), F(3, 5), F(2)),
            (F(3, 5), F(1), F(-2)),
        ),
        (
            (F(0), F(1, 7), F(5)),
            (F(1, 7), F(4, 7), F(-2)),
            (F(4, 7), F(6, 7), F(4)),
            (F(6, 7), F(1), F(1)),
        ),
    )

    covariance_checks = 0
    largest_ratio = F(0)
    for word in words:
        for function in functions:
            for exponent in range(5):
                frequency = 13**exponent
                error, bound = covariance_bound(function, word, frequency)
                require(error <= bound, "two-BV covariance bound")
                if bound:
                    largest_ratio = max(largest_ratio, error / bound)
                covariance_checks += 1

    word = words[0]
    alpha = word_mean(word)
    drift = functions[1]
    service = (
        (F(0), F(1, 3), F(2)),
        (F(1, 3), F(1), F(1)),
    )
    drift_total = mean_step(drift)
    service_total = mean_step(service)
    require(drift_total == 1 and service_total == F(4, 3), "control means")
    threshold_drift = (
        word_variation(word) * variation_step(drift)
        / (6 * alpha * abs(drift_total))
    )
    threshold_service = (
        word_variation(word) * variation_step(service)
        / (6 * alpha * service_total)
    )
    require(threshold_drift == 8, "drift threshold")
    require(threshold_service == F(3, 2), "service threshold")

    frequency = 13
    retained_drift = integrate_pulled_word(drift, word, frequency)
    retained_service = integrate_pulled_word(service, word, frequency)
    require(retained_drift >= alpha * abs(drift_total) / 2, "retained drift")
    require(retained_service >= alpha * service_total / 2, "retained service")

    require(4 * 63001 == 252004, "THM-2459 drift invoice")
    require(2 * 16384 == 32768, "THM-2459 service invoice")
    require(32768 * 2028 == 66453504, "worded root-coefficient invoice")

    half_word = ((F(0), F(1, 2)),)
    small_clock_drift = (
        (F(0), F(1, 2), F(0)),
        (F(1, 2), F(1), F(1)),
    )
    require(mean_step(small_clock_drift) == F(1, 2), "small-clock hostile mean")
    require(
        integrate_pulled_word(small_clock_drift, half_word, 1) == 0,
        "small-clock hostile",
    )

    zero_mean_drift = (
        (F(0), F(1, 2), F(1)),
        (F(1, 2), F(1), F(-1)),
    )
    whole_word = ((F(0), F(1)),)
    require(mean_step(zero_mean_drift) == 0, "zero-mean hostile mean")
    for exponent in range(5):
        require(
            integrate_pulled_word(zero_mean_drift, whole_word, 13**exponent) == 0,
            "zero-mean hostile survives",
        )

    print("THM-2466 DELAYED-WORD DRIFT/SERVICE AUDIT")
    print(f"covariance_checks={covariance_checks};largest_error_to_bound={largest_ratio}")
    print(f"word_density={alpha};word_variation={word_variation(word)}")
    print(
        "control_totals="
        f"drift:{drift_total},service:{service_total};"
        f"thresholds:{threshold_drift},{threshold_service}"
    )
    print(
        "frequency_13_retention="
        f"drift:{retained_drift},service:{retained_service}"
    )
    print(
        "thm2459_worded_floors="
        "alpha^2*D0/252004,alpha*M0/32768,root:alpha*M0/66453504"
    )
    print("small_clock_hostile=PASS;zero_mean_hostile=PASS")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
