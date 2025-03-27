from setuptools import setup, find_packages

setup(
    name='vachoppy',
    version='1.0.0',
    description='Python package for analyzing vacancy hopping mechanism',
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
        'mpi4py',
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
)
