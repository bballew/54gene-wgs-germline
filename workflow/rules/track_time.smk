rule track_start_time:
    output:
        "results/tat_tracking/start_time.txt",
    priority: 10000
    shell:
        "date +%s > {output}"
