from setuptools import setup

setup(
    name='KEGGX',
    packages=['KEGGX'],
    package_dir={'KEGGX': 'KEGGX'},
    package_data={'KEGGX': ['KEGG_compound_ids.txt']},
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
        "matplotlib==2.0.0", 
        "seaborn", 
        "requests"
    ]
)