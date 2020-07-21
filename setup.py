from setuptools import setup

from cellpose_filter.cellpose_filter import description

setup(
    name = "CellPose filter",
    packages = ["cellpose_filter"],
    version = 0.2,
    description = description,
    long_description = description,
    long_description_content_type='text/markdown',
    url='https://git.ist.ac.at/csommer/cellpose_filter',
    entry_points = {'console_scripts': [
            'cellpose_filter=cellpose_filter.cellpose_filter:main',
        ]},
    author = "Christoph Sommer",
    author_email = "christoph.sommer23@gmail.com",
    install_requires=[
            "scikit_image",
            "tifffile",
            "pandas",
            "matplotlib",
            "mpldatacursor",
            "openpyxl"
                    ]
    )