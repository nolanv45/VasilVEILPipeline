import sys

def parse_contig_orf(query_id):
    parts = query_id.split('_')
    contig_id = '_'.join(parts[:-3])
    return contig_id, query_id

def get_signature(target_id):
    return target_id.split('_')[0]

def process_files(file1_path, file2_path):
    primary_mappings = {}
    results = []
    unique_file2_orfs = set()
    
    # Process first file and store all its Query_IDs
    file1_queries = set()
    with open(file1_path, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            query_id = parts[0]
            target_id = parts[1]
            
            file1_queries.add(query_id)
            contig_id, orf_id = parse_contig_orf(query_id)
            signature = get_signature(target_id)
            
            primary_mappings[query_id] = signature
            # Add "Initial" for entries from first file
            results.append((contig_id, orf_id, signature, "Initial"))
    
    # Process second file
    with open(file2_path, 'r') as f:
        next(f)
        for line in f:
            parts = line.strip().split('\t')
            query_id = parts[0]
            target_id = parts[1]
            
            # If this Query_ID isn't in file1, add it to unique ORFs
            if query_id not in file1_queries:
                unique_file2_orfs.add(query_id)
            
            # Skip if query_id already exists in primary mappings
            if query_id in primary_mappings:
                continue
                
            # Look up the target_id in primary mappings
            if target_id in primary_mappings:
                contig_id, orf_id = parse_contig_orf(query_id)
                signature = primary_mappings[target_id]
                # Add "Recursive" for entries from second file
                results.append((contig_id, orf_id, signature, "Recursive"))
    
    return results, unique_file2_orfs

def main():
    if len(sys.argv) != 4:
        print("Usage: python script.py <top_hit_evalue_initial_search.tsv> <top_hit_evalue_recursive_search.tsv> <output_file>")
        sys.exit(1)
        
    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_file = sys.argv[3]
    
    # Process files
    results, unique_file2_orfs = process_files(file1_path, file2_path)
    
    # Write results to file
    with open(output_file, 'w') as f:
        # Updated header with new column
        f.write("Genome_ID\tORF_ID\tIdentified\tsignature\n")
        
        # Write data with identification source
        for contig_id, orf_id, signature, identified in results:
            f.write(f"{contig_id}\t{orf_id}\t{identified}\t{signature}\n")
    
    print(f"\nResults have been saved to: {output_file}")
    
    # Print unique ORFs from file 2
#    print("\nUnique ORFs found in file 2 (not present in file 1):")
#    for orf in sorted(unique_file2_orfs):
#        print(orf)
#    print(f"\nTotal number of unique ORFs in file 2: {len(unique_file2_orfs)}")

if __name__ == "__main__":
    main()
