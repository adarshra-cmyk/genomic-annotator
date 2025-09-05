from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open("requirements.txt", "r") as fh:
    requirements = [line.strip() for line in fh if line.strip()]

setup(
    name="genomic-annotator",
    version="1.0.0",
    author="Adarsh Ramamurthy",
    author_email="ramamurthyadarsh8@gmail.com",
    description="Genetic Annotator SDK",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/adarshra-cmyk/genomic-annotator",
    packages=find_packages(),
    python_requires=">=3.7",
    install_requires=requirements,
    include_package_data=True,
)
