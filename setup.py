import os

from setuptools import find_packages, setup

INSTALL_REQUIRES = ["numpy", "scikit-learn"]
EXTRAS_REQUIRE = {"cli":["click", "pandas"]}


def read_version():
    version_py = os.path.join(os.path.dirname(__file__), "python", "version.py")
    with open(version_py) as handle:
        version = handle.read().strip().split("=")[-1].replace('"', "").strip()
    return version


def read_readme():
    with open("README.md") as handle:
        readme = handle.read()
    return readme


setup(
    name="knn_smooth",
    version=read_version(),
    packages=["knn_smooth"],
    package_dir={"knn_smooth":"python"},
    license="MIT",
    install_requires=INSTALL_REQUIRES,
    extras_require=EXTRAS_REQUIRE,
    python_requires=">=3.6",
    entry_points={"console_scripts": ["knn_smooth=knn_smooth.cli:main"],},
    long_description=read_readme(),
    long_description_content_type="text/markdown",
)
