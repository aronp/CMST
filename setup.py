from setuptools import setup, find_packages

# Read the README for the long description on PyPI
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="cmst-window",
    version="0.1.0",
    author="Aron Palmer",
    author_email="palmer.aron+cmst@gmail.com", 
    description="The analytically stable, zero-preserving, hyper-flat window.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/aronp/CMST",
    project_urls={
        "Bug Tracker": "https://github.com/aronp/CMST/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License ::OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Development Status :: 4 - Beta",
    ],
    # This section is critical for the 'src' layout
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.18.0",
    ],
)
