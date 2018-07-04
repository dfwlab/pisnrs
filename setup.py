import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pisnrs",
    version="0.0.2",
    author="Dingfeng Wu",
    author_email="dfw_bioinfo@126.com",
    description="A package for predict scaffold of NR inhibitors",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pypa/pisnrs",
    packages=setuptools.find_packages(),
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ),
    install_requires=[
        'scikit-learn>=0.19.1',
        'pandas>=0.23.1'
    ],
)