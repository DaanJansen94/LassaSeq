from setuptools import setup, find_packages

setup(
    name="lassaseq",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "biopython>=1.81",
    ],
    entry_points={
        'console_scripts': [
            'lassaseq=lassaseq.lassaseq:cli_main',
        ],
    },
    author="Daan Jansen",
    description="A tool for filtering and analyzing Lassa virus sequences",
) 