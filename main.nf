nextflow.enable.dsl=2

process sayHello {
    output:
    stdout

    script:
    """
    echo "Hello, Nextflow!"
    """
}

workflow {
    // Run the process and capture its output
    ch = sayHello()
    
    // View the output by printing to the console
    ch.view()
}
