from setuptools import setup

setup(
    name='KEGGX',
    packages=['KEGGX'],
    package_dir={'KEGGX': 'src'},
    include_package_data=True,
    version='0.1.0',
    url='https://github.com/iamjli/KEGGX',
    python_requires='>=3.5',
    classifiers=[
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    license='MIT',
    author='iamjli',
    author_email='iamjli@mit.edu',
    description='',
    install_requires=[
        "pandas>=0.20.1", 
        "numpy>=1.12.0",
        "networkx>=2.0",
    ]
)