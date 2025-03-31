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
    sayHello.view()
}
