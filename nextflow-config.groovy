// Nextflow configuration file

// Process configuration
process {
    // Default process executor
    executor = 'local'
    
    // Default number of CPUs
    cpus = params.threads
    
    // Error handling - retry up to 3 times
    errorStrategy = { task.exitStatus in [1,143,137,104,134,139] ? 'retry' : 'finish' }
    maxRetries = 3
    maxErrors = '-1'
}

// Environment configuration
env {
    // Set locale to handle special characters
    LC_ALL = "en_US.UTF-8"
}

// Executor configuration (customize based on your infrastructure)
executor {
    $local {
        cpus = 8
        memory = '16 GB'
    }
    
    // Uncomment and customize if using SLURM cluster
    /*
    $slurm {
        queueSize = 100
        jobName = { "rnaseq_${task.index}" }
        queue = 'standard'
    }
    */
}

// Enable container usage (Docker or Singularity)
// Uncomment and customize for your environment
/*
process {
    container = 'your-container-image:latest'
}

singularity {
    enabled = true
    autoMounts = true
}
*/

// Report configuration
report {
    enabled = true
    file = "${params.output_dir}/reports/execution_report.html"
}

timeline {
    enabled = true
    file = "${params.output_dir}/reports/timeline.html"
}

trace {
    enabled = true
    file = "${params.output_dir}/reports/trace.txt"
}

// Profiles configuration
profiles {
    standard {
        process.executor = 'local'
    }
    
    cluster {
        process.executor = 'slurm'
        process.queue = 'standard'
    }
    
    cloud {
        process.executor = 'awsbatch'
        aws.region = 'us-east-1'
        aws.batch.cliPath = '/usr/local/bin/aws'
    }
}

// Manifest
manifest {
    name = 'RNA-seq Pipeline with Minimap2 and Oarfish'
    author = 'Your Name'
    description = 'Nextflow pipeline for RNA-seq analysis using Minimap2 and Oarfish'
    version = '1.0.0'
    mainScript = 'main.nf'
}
