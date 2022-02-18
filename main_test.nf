nextflow.enable.dsl = 2

Channel.
      from("a", "b", 3).
      collect().
      set{ test }

log.info test
