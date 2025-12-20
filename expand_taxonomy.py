import argparse
import csv
import sys
from Bio import Entrez

def get_descendants(taxid):
    """
    Fetches all descendant TaxIDs for a given TaxID using NCBI Entrez.
    """
    try:
        # Search for all items in the subtree of the given taxid
        handle = Entrez.esearch(db="taxonomy", term=f"txid{taxid}[Subtree]", RetMax=100000)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Error fetching descendants for TaxID {taxid}: {e}", file=sys.stderr)
        return []

def get_taxonomy_details(taxids):
    """
    Fetches ScientificName and Rank for a list of TaxIDs.
    Batches requests to avoid URI length issues.
    """
    details = {}
    batch_size = 200

    for i in range(0, len(taxids), batch_size):
        batch = taxids[i:i+batch_size]
        try:
            handle = Entrez.efetch(db="taxonomy", id=",".join(batch))
            records = Entrez.read(handle)
            handle.close()

            for record in records:
                t_id = record['TaxId']
                name = record.get('ScientificName', 'Unknown')
                rank = record.get('Rank', 'no rank')
                details[t_id] = {'ScientificName': name, 'Rank': rank}
        except Exception as e:
            print(f"Error fetching details for batch starting at index {i}: {e}", file=sys.stderr)

    return details

def main():
    parser = argparse.ArgumentParser(description="Expand taxonomy in CSV file.")
    parser.add_argument("--input_file", help="Input CSV file", default="human_pathogen_liverpool.csv")
    parser.add_argument("--output_file", help="Output CSV file", required=True)
    parser.add_argument("--email", help="Email for NCBI Entrez", default="jochen_schaefer@gmx.de")
    parser.add_argument("--api_key", help="NCBI Entrez API Key", default=None)

    args = parser.parse_args()

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    existing_taxids = set()

    rows = []

    def read_csv_with_fallback(filepath):
        encodings = ['utf-8', 'latin-1', 'cp1252']
        for enc in encodings:
            try:
                with open(filepath, 'r', newline='', encoding=enc) as f:
                    # Read entire file to force decoding check
                    return list(csv.reader(f, delimiter=';'))
            except UnicodeDecodeError:
                continue
            except Exception as e:
                raise e
        raise ValueError(f"Could not decode file {filepath} with encodings: {encodings}")

    try:
        all_rows = read_csv_with_fallback(args.input_file)
        if not all_rows:
            print("File is empty.")
            sys.exit(1)

        header = all_rows[0]

        # Find indices
        try:
            taxid_idx = header.index("TaxId")
            organism_idx = header.index("Organism")

            # Verify indices based on name might be ambiguous for "Rank".
            # Let's find all indices for "Rank"
            rank_indices = [i for i, x in enumerate(header) if x == "Rank"]
            if len(rank_indices) >= 2:
                rank_numeric_idx = rank_indices[0]
                rank_tax_idx = rank_indices[1]
            else:
                # Fallback or error
                rank_numeric_idx = 0
                rank_tax_idx = 4 # Guessing from example if header names don't match exactly

        except ValueError as e:
            print(f"Error parsing header: {e}")
            sys.exit(1)

        for row in all_rows[1:]:
            rows.append(row)
            if len(row) > taxid_idx:
                existing_taxids.add(row[taxid_idx])

    except FileNotFoundError:
        print(f"File {args.input_file} not found.")
        sys.exit(1)
    except ValueError as e:
        print(e)
        sys.exit(1)

    # Process
    new_rows = []

    for row in rows:
        parent_taxid = row[taxid_idx]
        print(f"Processing TaxID: {parent_taxid}")

        descendants = get_descendants(parent_taxid)

        # Filter out already existing
        to_fetch = [d for d in descendants if d not in existing_taxids]

        if not to_fetch:
            continue

        print(f"Found {len(to_fetch)} new descendants for {parent_taxid}")

        # Fetch details for new descendants
        details_map = get_taxonomy_details(to_fetch)

        for d_taxid in to_fetch:
            if d_taxid not in details_map:
                continue

            info = details_map[d_taxid]
            new_row = list(row) # Copy parent row

            # Update fields
            # TaxId
            new_row[taxid_idx] = d_taxid
            # Organism
            new_row[organism_idx] = info['ScientificName']
            # Rank (Taxonomy)
            new_row[rank_tax_idx] = info['Rank']

            new_rows.append(new_row)
            existing_taxids.add(d_taxid) # Add to set to avoid duplicates if processed again

    # Write output
    with open(args.output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter=';')
        writer.writerow(header)
        writer.writerows(rows)
        writer.writerows(new_rows)

    print(f"Finished. Wrote {len(rows) + len(new_rows)} rows to {args.output_file}")

if __name__ == "__main__":
    main()
