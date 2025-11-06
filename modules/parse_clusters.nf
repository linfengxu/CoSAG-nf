#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PARSE_FILTERED_CLUSTERS {
    container "${params.python.container}"
    input:
    path filtered_json
    
    output:
    path "cluster_*.txt", emit: cluster_files
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    
    # 读取过滤后的JSON文件
    with open('${filtered_json}', 'r') as f:
        data = json.load(f)
    
    # 为每个cluster创建一个文件
    for cluster in data['clusters']:
        cluster_id = cluster['cluster_id']
        filename = f"cluster_{cluster_id}.txt"
        
        with open(filename, 'w') as f:
            f.write(f"cluster_id\\t{cluster_id}\\n")
            f.write(f"cluster_size\\t{cluster['cluster_size']}\\n")
            f.write("members\\n")
            
            for member in cluster['members']:
                f.write(f"{member['sag_id']}\\t{member['read1']}\\t{member['read2']}\\n")
    """
}

process PARSE_SUBCLUSTER_FILES {
    container "${params.python.container}"
    input:
    path subcluster_files
    
    output:
    path "subcluster_*.txt", emit: subcluster_files
    
    script:
    """
    #!/usr/bin/env python3
    
    import json
    import glob
    
    # 处理所有subcluster JSON文件
    json_files = glob.glob('*.json')
    
    for json_file in json_files:
        with open(json_file, 'r') as f:
            data = json.load(f)
        
        original_cluster_id = data['original_cluster_id']
        
        for subcluster in data['subclusters']:
            subcluster_id = subcluster['subcluster_id']
            filename = f"subcluster_{subcluster_id}.txt"
            
            with open(filename, 'w') as f:
                f.write(f"original_cluster_id\\t{original_cluster_id}\\n")
                f.write(f"subcluster_id\\t{subcluster_id}\\n")
                f.write("members\\n")
                
                for member in subcluster['members']:
                    f.write(f"{member['sag_id']}\\t{member['read1']}\\t{member['read2']}\\n")
    """
}