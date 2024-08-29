import setuptools

# Get some values from the setup.cfg
try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

conf = ConfigParser()
conf.read(['setup.config'])
metadata = dict(conf.items("metadata"))

NAME = metadata['name']
VERSION = metadata['version']
RELEASE = 'dev' not in VERSION
AUTHOR = metadata["author"]
AUTHOR_EMAIL = metadata["author_email"]
LICENSE = metadata["license"]
DESCRIPTION = metadata["description"]

entry_points = {
    'console_scripts': [
        "kcrm_create_crmsk = scripts.kcrm_create_crmsk:main",
        "kcrm_group_famres = scripts.kcrm_group_frames:main",
        "kcwi_collapse = scripts.kcwi_collapse:main",
        "kcwi_combinestd = scripts.kcwi_combinestd:main",
        "kcwi_flatten_cube = scripts.kcwi_flatten_cube:main",
        "kcwi_gen_log = scripts.kcwi_gen_log:main",
        "kcwi_makemask_medfilter = scripts.kcwi_makemask_medfilter:main",
        "kcwi_masksky_ds9_2d = scripts.kcwi_masksky_ds9_2d:main",
        "kcwi_masksky_ds9_thum = scripts.kcwi_masksky_ds9_thum:main",
        "kcwi_masksky_ds9 = scripts.kcwi_masksky_ds9:main",
        "kcwi_medfilter = scripts.kcwi_medfilter:main"
    ],
    'gui_scripts': [
        "kcwi_viewer = gui.kcwi_viewer:main"
    ]
}

setuptools.setup(name=NAME,
      provides=NAME,
      version=VERSION,
      license=LICENSE,
      description=DESCRIPTION,
      long_description=open('README.md').read(),
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      packages=setuptools.find_packages(),
      package_data={'': ['data/extin/*', 'data/stds/*']},
      entry_points=entry_points,
      install_requires=[
        'astropy==4.3.1',
        'argparse',
        'matplotlib==3.5.0',
        'numpy==1.21.2',
        'scipy==1.7.1']
)


