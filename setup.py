from setuptools import setup, find_packages

setup(
    name='BioHelperFunc',
    version='0.1',
    packages=find_packages(),
    author="JustABiologist",
    author_email="stab.me.papi@gmail.com",
    description="Short functions that might help out someone idk",
    url="https://github.com/JustABiologist/BioHelperFunctions",
    install_requires=[
        "biopython",
        "pandas",
        "numpy",
        "matplotlib",
        "tkinter"
    ],
)