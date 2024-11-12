from setuptools import setup, find_packages

setup(
    name='vachoppy',
    version='0.0.2',
    packages = find_packages(),
    entry_points={
        'console_scripts': [
            'vachoppy=vachoppy.main:main',
        ],
    },
)
