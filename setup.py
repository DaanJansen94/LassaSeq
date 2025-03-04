from setuptools import setup, find_packages

setup(
    name="lassaseq",
    version="0.1.0",
    packages=find_packages(),
    package_data={
        'lassaseq': ['lineages/*'],
    },
    include_package_data=True,
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
    description="Tool for downloading and analyzing Lassa virus sequences",
    author="Daan Jansen",
    author_email="jansendaan94@gmail.com",
    url="https://github.com/DaanJansen94/LassaSeq"
)
