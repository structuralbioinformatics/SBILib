from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()


setup(
    name='SBILib',
    version='0.3.3',
    license='MIT',
    author="Patrick Gohl",
    author_email='patrick.gohl@upf.edu',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    include_package_data=True,
    long_description=readme,
    long_description_content_type="text/markdown",
    url='https://github.com/structuralbioinformatics/SBI',
    keywords='Structural Bioinformatics, Loops',
    install_requires=[
          'numpy','pandas','pynion','Requests','scipy','six'
      ],

)




