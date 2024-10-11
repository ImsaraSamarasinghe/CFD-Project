from setuptools import setup, find_packages

setup(
    name='cfd_simulation',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'matplotlib',
    ],
    entry_points={
        'console_scripts': [
            'run_simulation=cfd_simulation.main:main',
        ],
    },
    description='A Python package for solving cavity flow problems using finite difference methods.',
    author='ImsaraSamarasinghe',
    author_email='imsara256@gmail.com',
    url='https://github.com/ImsaraSamarasinghe?tab=repositories',
    license='MIT',
)
