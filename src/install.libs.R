# ref : h/t @Zilong-li https://github.com/Zilong-Li/vcfppR/blob/main/src/install.libs.R
## Install shared objects and static libraries
files <- Sys.glob(paste0("*", SHLIB_EXT))

## Also copy the bcftools static library
files <- c(files, dir("bcftools-1.21/", pattern = "a$", full.names = TRUE))

## Create destination directory
dest <- file.path(R_PACKAGE_DIR, paste0("libs", R_ARCH))
dir.create(dest, recursive = TRUE, showWarnings = FALSE)

## Copy the shared objects and static libraries
file.copy(files, dest, overwrite = TRUE)

## Copy header files
headers <- c(
    list.files("bcftools-1.21/", pattern = "\\.h$", full.names = TRUE)
)
headers_dest <- file.path(R_PACKAGE_DIR, "include/bcftools")
dir.create(headers_dest, recursive = TRUE, showWarnings = FALSE)
file.copy(headers, headers_dest, overwrite = TRUE)

## Print message for debugging purposes
cat("Installed libraries:", paste(files, collapse = ", "), "\n")
cat("Installed headers:", paste(headers, collapse = ", "), "\n")
