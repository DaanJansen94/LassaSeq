from setuptools import setup, find_packages

setup(
    name="lassaseq",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'biopython',
        'requests',
        'tqdm',
    ],
    entry_points={
        'console_scripts': [
            'lassaseq=lassaseq.lassaseq:cli_main',
        ],
    },
    description="Tool for downloading Lassa virus sequences",
    author="Your Name",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/lassaseq"
)
