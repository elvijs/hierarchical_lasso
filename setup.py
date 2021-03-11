"""
Hierarchical Lasso models.

This package implements Hierarchical Lasso regression models as discussed in
https://arxiv.org/abs/2001.07778.
"""

from distutils.core import setup

REPO_URL = "https://github.com/elvijs/hierarchical_lasso"


def setup_package():
    setup(
        name="Hierarchical Lasso",
        version="0.0.1",
        description="Hierarchical Lasso algorithm.",
        long_description=__doc__,
        author="Hugo Maruri-Aguilar",
        author_email="h.maruri-aguilar@qmul.ac.uk",
        maintainer="Elvijs Sarkans",
        maintainer_email="Elvijs.Sarkans@gmail.com",
        url=REPO_URL,
        download_url=REPO_URL,
        license="MIT",
        packages=['src/hierarchical_lasso'],
        package_data={},
        install_requires=[
            "scipy>=1.6",
        ],
        python_requires=">=3.5",
    )


if __name__ == '__main__':
    setup_package()
