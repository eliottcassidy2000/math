"""
Setup file for the cayley-delannoy package.

Install: pip install -e .
Usage:   from cayley_delannoy import Q, gk, cv2, W
CLI:     python -m cayley_delannoy.tournament_test --n 7 --h 189
"""

from setuptools import setup

setup(
    name="cayley-delannoy",
    version="0.1.0",
    description="Cayley-Delannoy tournament theory: ranking significance via lattice path combinatorics",
    author="Eliott Cassidy",
    py_modules=["cayley_delannoy", "tournament_test", "ab_test_ranker"],
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "tournament-test=tournament_test:main",
            "ab-test-ranker=ab_test_ranker:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
)
