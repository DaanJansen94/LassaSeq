from setuptools import setup, find_packages

setup(
    name="lassaseq",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'biopython>=1.80',
        'numpy>=1.20.0',
        'requests>=2.25.0',
        'tqdm>=4.50.0',
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
