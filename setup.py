from setuptools import setup, find_packages

setup(
    name="dardar_tool",
    version="1.0.0",
    description="Tool to download, extract and plot dardar data.",
    author="Arttu Väisänen",
    packages=find_packages(),
    package_data={"dardar_tool": ["source/*.nc"]},
    python_requires="==3.10.18",
    install_requires=[
        "numpy==2.2.6",
        "matplotlib==3.10.3",
        "netCDF4==1.7.2",
        "cartopy==0.24.0",
        "cryptography==45.0.3",
        "paramiko==3.5.1",
        "cftime==1.6.4",
        ],
    )
