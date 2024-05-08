## Species Distribution Modeling (SDM) Code

This R code repository contains scripts for conducting Species Distribution Modeling (SDM) analyses. The code is developed as part of a master thesis project and aims to predict species distribution in a glacier foreland through SDM and to assess its transferability between two adjacent glacier forelands with comparable environmental conditions.


## Installation

The dependencies required to run the code and provided in the "SDM.R" file, in the "Load libraries" section. 

You'll need to download an older version of the `biomod2` package, as the code relies on specific functionality provided by this version. Here's how to proceed:

1. Download the package archive for version 4.2-4 of the `biomod2` package from the following URL: [biomod2_4.2-4.tar.gz](https://cran.r-project.org/src/contrib/Archive/biomod2/biomod2_4.2-4.tar.gz).

2. Save the downloaded file to your local computer, in the repository directory.

3. Use the `install.packages()` function in R to install the package from your local directory. Here's an example:

```R
# Replace "C:/path/to/biomod2_4.2-4.tar.gz" with the path where you saved the downloaded file.
install.packages("C:/path/to/biomod2_4.2-4.tar.gz", repos = NULL, type = "source")

# Load the installed package
library(biomod2)


Note: If you encounter an error mentioning 00LOCK-biomod2, it means that there is a directory lock preventing R from modifying the library directory. To resolve this issue:

1. Open File Explorer and navigate to your R library directory, typically located at C:/Users/.../AppData/Local/R/win-library/4.3.

2. Look for a folder named 00LOCK-biomod2.

3. Delete this folder or move the lock file outside of the R library directory.

4. After removing the lock, retry installing the package using the command:

install.packages("C:/path/to/biomod2_4.2-4.tar.gz", repos = NULL, type = "source")

If you continue to encounter issues, ensure that no other processes are using the R library directory. Close any other programs that might be accessing R libraries simultaneously. If the problem persists, consider restarting your computer before attempting to install the package again.
