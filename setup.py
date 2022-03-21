import setuptools

setuptools.setup(
    name='deltapersistence',
    version='1.0',
    author='Joshua Mirth, Johnathan Bush',
    author_email='joshua.mirth@gmail.com',
    description='Persistent homology tools for NSF-DELTA.',
    url='https://gitlab.com/jrmirth/deltapersistence',
    project_urls={
        "Documentation": "https://jrmirth.gitlab.io/deltapersistence/",
        "Delta Homepage": "https://delta-topology.org/",
    },
    license='GNU GPLv3',
    packages=setuptools.find_packages(),
    scripts=['persistence_script.py'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPL License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=['numpy'],
    zip_safe=False
)
