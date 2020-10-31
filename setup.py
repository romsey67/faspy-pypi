import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="faspy-pkg-ROSLI-MOHD-SANI", # Replace with your own username
    version="0.0.2",
    author="ROSLI MOHD SANI",
    author_email="romsey67@gmail.com",
    description="A python package for financial instruments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/romsey67/faspy.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
