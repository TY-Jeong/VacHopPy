from setuptools import setup, find_packages

__version__="0.0.0"

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

setup(
    name='vachoppy',
    version=__version__,
    description='Python package for analyzing vacancy hopping mechanism',
    long_description = long_description,
    long_description_content_type='text/markdown',
    author='TY-Jeong',
    author_email='helianthus312@gmail.com',
    url='https://github.com/TY-Jeong/VacHopPy',
    packages = find_packages(),
    install_requires=[
        'numpy>=1.26.4',
        'tqdm>=4.67.1',
        'colorama',
        'matplotlib>=3.10.0',
        'scipy',
        'tabulate',
        'pymatgen>=2024.6.10'
    ],
    extras_require={
        'parallel': ['mpi4py']
    },
    python_requires='>=3.10',
    keywords=['vachoppy', 'vacancy', 'hopping'],
    entry_points={
        'console_scripts': [
            'vachoppy=vachoppy.main:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics"
    ],
)
