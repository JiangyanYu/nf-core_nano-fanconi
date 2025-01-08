process TABIX_TABIX_COPY {
    input:
    path tbi_file // Input file to copy
    path target_dir // Directory to copy into

    output:
    path "${target_dir}/${tbi_file.name}" // The copied file in the target directory

    script:
    """
    cp $tbi_file $target_dir/
    """
}
