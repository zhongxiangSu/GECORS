# coding utf8
import setuptools
from gecors.versions import get_versions

with open('README.md') as f:
    LONG_DESCRIPTION = f.read()

setuptools.setup(
    name="GECORS",
    version=get_versions(),
    author="zhongxiang Su",
    author_email="suzhongxiang@mail.kib.ac.cn",
    description="Get diploid genome cricle map",
    long_description=LONG_DESCRIPTION,
    long_description_content_type='text/markdown',
    url="https://github.com/zhongxiangSu/GECORS",

    entry_points={
        "console_scripts": ["GECORS = gecors.cli:main"]
    },    

    packages=setuptools.find_packages(),

    install_requires=[
        "matplotlib>=3.0.3",
    ],

    python_requires='>=3.5',
)