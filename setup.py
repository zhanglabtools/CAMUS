from setuptools import Command, find_packages, setup


setup(
    name = "CAMUS",
    version = "1.0.0",
    description = "Cross-datasets annotation modeling with universal reference data and method selection",
    url = "https://github.com/Zhanglabtools/CAMUS",
    author = "Qunlun Shen",
    author_email = "knotnet@foxmail.com",
    license = "MIT",
    package_dir={'': 'CAMUS'},
    install_requires = ["requests",],
    zip_safe = False,
    include_package_data = True,
)