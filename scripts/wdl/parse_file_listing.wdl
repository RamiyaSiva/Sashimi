## This WDL script accepts a string with a list of file paths and parses the list
##     to return and Array of files
##
## Requirements/expectations :
## - String containing a list of file paths in the format:
##     ["file1.abc", "file2.xyz", ...]

# WORKFLOW DEFINITION
workflow ParseFileListing {
    String source_files_str

    # Parse file list string and return list of files
    call ParseFileList {
        input: file_list=source_files_str
    }

    output {
        Array[String] source_files = ParseFileList.files
    }
}

# Parse file list string and return list of files
task ParseFileList {
    String file_list

    command {
        # Convert delimited string to list of files
        echo "${file_list}" | tr -d '[] "' | sed 's/,/\n/g'
   }

    runtime {
        docker: 'ubuntu:latest'
    }
  
    output {
        Array[String] files = read_lines(stdout())
    }
}
