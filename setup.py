from setuptools import Extension, setup
module=Extension("symnmfsp", sources=["symnmfmodule.c", "symnmf.c"])
setup(name="symnmfsp",
      version="1.0",
      description="Python wrapper for custom C extension",
      ext_modules=[module]
      )
