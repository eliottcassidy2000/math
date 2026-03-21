from setuptools import setup, find_packages

setup(
    name="tournament-toolkit",
    version="0.1.0",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=["numpy>=1.20"],
    author="Tournament Theory Research Project",
    description="Mathematical tools for ranking, fraud detection, and AI diagnostics based on tournament theory",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/eliottcassidy2000/math",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Mathematics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
    ],
)
