with open('README.md') as readme_file:
    readme = readme_file.read()

exec(open('bgc_mcl/version.py').read()) # loads __version__

setup(name='bgc_mcl',
      version=__version__,
      author='fernando garcia-guevara',
      description='Cluster BGC by gene content',
      long_description=readme,
      license='GPL3+',
      keywords="Biosynthetic Gene Cluster MCL",
      packages=find_packages(exclude='docs'),
      #install_requires=('biopython >=1.64'),
      url='https://github.com/garciafertson/BGC_MCL',
      scripts=['bin/BGCfind_hmm.py'],
)
                                    
