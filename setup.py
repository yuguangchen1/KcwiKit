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
        "kcrm_create_crmsk = kcwikit.scripts.kcrm_create_crmsk:main",
        "kcrm_group_frames = kcwikit.scripts.kcrm_group_frames:main",
        "kcwi_collapse = kcwikit.scripts.kcwi_collapse:main",
        "kcwi_combinestd = kcwikit.scripts.kcwi_combinestd:main",
        "kcwi_flatten_cube = kcwikit.scripts.kcwi_flatten_cube:main",
        "kcwi_gen_log = kcwikit.scripts.kcwi_gen_log:main",
        "kcwi_makemask_medfilter = kcwikit.scripts.kcwi_makemask_medfilter:main",
        "kcwi_makemask = kcwikit.scripts.kcwi_makemask:main",
        "kcwi_masksky_ds9_2d = kcwikit.scripts.kcwi_masksky_ds9_2d:main",
        "kcwi_masksky_ds9_thum = kcwikit.scripts.kcwi_masksky_ds9_thum:main",
        "kcwi_masksky_ds9_zap = kcwikit.scripts.kcwi_masksky_ds9_zap:main",
        "kcwi_masksky_ds9 = kcwikit.scripts.kcwi_masksky_ds9:main",
        "kcwi_medfilter = kcwikit.scripts.kcwi_medfilter:main",
        "kcwi_masksky_ds9_bcrr = kcwikit.scripts.kcwi_masksky_ds9_bcrr:main",
        "kcwi_gen_skyfile = kcwikit.scripts.kcwi_gen_skyfile:main"
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
      package_data={'': ['data/extin/*', 'data/stds/*']}, # These should belong to KSkyWizard, but are tmp kept here for safety purposes
      entry_points=entry_points,
      install_requires=[
        'astropy',
        'argparse',
        'matplotlib',
        'numpy',
        'scipy',
        'regions',
        'fpdf',
        'bokeh',
        'colorcet',
        'tqdm',
        'MontagePy']
)



