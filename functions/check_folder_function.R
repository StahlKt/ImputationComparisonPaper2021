# check directories and create them, if they don't exist yet 

check_and_create_folders <- function(directory) {
    if (!dir.exists(directory)) {
        print("folder created")
        dir.create(directory, recursive = TRUE)
        return()
    } else {
        print("folder exists")
        return()
    }
}

