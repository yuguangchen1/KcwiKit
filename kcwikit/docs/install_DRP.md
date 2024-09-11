# Installing the modified DRP

**Note**: A new DRP version that supports `Python >= 3.11` is currently being prepared by the observatory. Some of the improvements on this page will be included in the new DRP. 

**Disclaimer**: Some features introduced in this modified version has not been sufficiently tested (which is why it is not yet merged into the official version). We are happy to receive bug reports, but we cannot guarantee that this version would work for your observations. So use with caution!

There are two methods to install the modified DRP. 

1. If you have installed the official DRP with `Python ~= 3.8`, you can download the modified files from the [pyDRP](../pyDRP/) directory in this repository and replace the older ones with them. Then, rerun the installation script (`python setup.py install`). 

2. Alternatively, we keep a separate version of DRP at `https://github.com/yuguangchen1/KCWI_DRP/tree/KCWIKit`. Download this updated version.

    ```bash
    git clone https://github.com/yuguangchen1/KCWI_DRP.git
    git checkout KCWIKit
    ```

    And install this updated version. 

    ```bash
    python setup.py install
    ```


