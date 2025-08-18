from setuptools import setup  # noqa: D100

with open("README.md", encoding="utf-8") as f:  # noqa: FURB101
    long_description = f.read()


setup(
    name="cobrak",
    description="Implementation of COBRA-k",
    long_description=open("README.md", encoding="utf8").read(),  # noqa: SIM115
    long_description_content_type="text/asciidoc",
    packages=["cobrak"],
    package_dir={"cobrak": "cobrak"},
    # package_data={'cobrak': [
    #     'data/*.svg',
    # ]}
)
