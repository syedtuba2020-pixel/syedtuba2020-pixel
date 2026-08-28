"""
MINI PROJECT
Biomedical Data Collector using Open Targets and UniProt APIs

Name: Tuba Durdana

OPTION 1 — Disease Target Data Collector
OPTION 2 — Protein Information Collector
"""

"""
OPTION 1 — Disease Target Data Collector
Biomedical Data Collector using the Open Targets Platform GraphQL API.

Input: disease name
Output: disease_targets.csv
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd
import requests

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

SEARCH_QUERY = """
query DiseaseSearch($queryString: String!) {
  search(queryString: $queryString, entityNames: ["disease"]) {
    hits {
      id
      name
    }
  }
}
"""

TARGET_QUERY = """
query AssociatedTargets($efoId: String!, $index: Int!, $size: Int!) {
  disease(efoId: $efoId) {
    id
    name
    associatedTargets(page: {index: $index, size: $size}) {
      count
      rows {
        target {
          id
          approvedSymbol
          approvedName
          biotype
        }
        score
      }
    }
  }
}
"""


def post_graphql(query, variables):
    response = requests.post(
        API_URL,
        json={"query": query, "variables": variables},
        timeout=30,
        headers={"Content-Type": "application/json"},
    )
    response.raise_for_status()
    payload = response.json()

    if payload.get("errors"):
        messages = "; ".join(
            error.get("message", "Unknown GraphQL error")
            for error in payload["errors"]
        )
        raise RuntimeError(messages)

    return payload.get("data", {})


def find_disease(disease_name):
    data = post_graphql(SEARCH_QUERY, {"queryString": disease_name})
    hits = data.get("search", {}).get("hits", [])

    if not hits:
        raise ValueError(
            f"No disease was found for '{disease_name}'. "
            "Try a more specific disease name or provide an EFO/MONDO ID."
        )

    # Prefer an exact name match when available.
    exact = [h for h in hits if h.get("name", "").lower() == disease_name.lower()]
    return (exact or hits)[0]


def get_associated_targets(efo_id, page_size=100, max_rows=500):
    rows = []
    index = 0
    disease_name = ""

    while len(rows) < max_rows:
        data = post_graphql(
            TARGET_QUERY,
            {"efoId": efo_id, "index": index, "size": min(page_size, max_rows - len(rows))},
        )
        disease = data.get("disease")
        if not disease:
            raise ValueError(f"Open Targets returned no disease for ID '{efo_id}'.")

        disease_name = disease.get("name", "")
        association = disease.get("associatedTargets") or {}
        page_rows = association.get("rows") or []
        if not page_rows:
            break

        for item in page_rows:
            target = item.get("target") or {}
            rows.append(
                {
                    "Gene Symbol": target.get("approvedSymbol"),
                    "Gene Name": target.get("approvedName"),
                    "Biotype": target.get("biotype"),
                    "Association Score": item.get("score"),
                    "Ensembl ID": target.get("id"),
                }
            )

        total = association.get("count", len(rows))
        if len(rows) >= total or len(page_rows) < page_size:
            break
        index += 1

    return disease_name, pd.DataFrame(rows)


def safe_filename(name):
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return name.strip("_") or "disease"


def collect_disease_targets(disease_name=None, disease_id=None, output="disease_targets.csv"):
    if disease_id:
        selected_id = disease_id.strip()
        selected_name = selected_id
    else:
        if not disease_name:
            raise ValueError("Disease name is required.")
        selected = find_disease(disease_name)
        selected_id = selected["id"]
        selected_name = selected["name"]

    resolved_name, df = get_associated_targets(selected_id)

    if df.empty:
        raise ValueError("No associated targets were returned.")

    # Highest association score first.
    df = df.sort_values(
        by="Association Score", ascending=False, na_position="last"
    ).reset_index(drop=True)

    df.to_csv(output, index=False)

    print(f"\nDisease: {resolved_name or selected_name}")
    print(f"Disease ID: {selected_id}")
    print(f"Targets retrieved: {len(df)}")
    print(f"Saved: {output}\n")
    print(df.head(10).to_string(index=False))
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve disease-associated targets from Open Targets."
    )
    parser.add_argument("disease", nargs="?", help="Disease name, e.g. Breast Cancer")
    parser.add_argument("--disease-id", help="Known EFO/MONDO disease ID")
    parser.add_argument(
        "-o", "--output", default="disease_targets.csv", help="CSV output path"
    )
    args = parser.parse_args()

    disease = args.disease
    if not disease and not args.disease_id:
        disease = input("Enter disease name (e.g. Breast Cancer): ").strip()

    if not disease and not args.disease_id:
        sys.exit("Please enter a disease name or provide --disease-id.")

    try:
        collect_disease_targets(disease, args.disease_id, args.output)
    except requests.RequestException as exc:
        sys.exit(f"Network/API error: {exc}")
    except Exception as exc:
        sys.exit(f"Error: {exc}")


if __name__ == "__main__":
    main()


"""
OPTION 2 — Protein Information Collector
Biomedical Data Collector using the UniProt REST API.

Input: protein/gene name
Output: <query>_proteins.csv
"""

import argparse
import re
import sys

import pandas as pd
import requests

API_URL = "https://rest.uniprot.org/uniprotkb/search"

FIELDS = "accession,protein_name,gene_names,organism_name,length"


def safe_filename(name):
    name = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
    return name.strip("_") or "protein"


def search_uniprot(query, size=10):
    params = {
        "query": f'"{query}"',
        "format": "json",
        "fields": FIELDS,
        "size": size,
    }

    response = requests.get(API_URL, params=params, timeout=30)
    response.raise_for_status()
    payload = response.json()

    records = []
    for item in payload.get("results", []):
        genes = item.get("genes") or []
        gene_names = []

        for gene in genes:
            primary = (gene.get("geneName") or {}).get("value")
            if primary:
                gene_names.append(primary)

            for synonym in gene.get("synonyms", []) or []:
                value = synonym.get("value")
                if value:
                    gene_names.append(value)

        records.append(
            {
                "Accession": item.get("primaryAccession"),
                "Protein Name": (item.get("proteinDescription") or {})
                .get("recommendedName", {})
                .get("fullName", {})
                .get("value"),
                "Gene": "; ".join(dict.fromkeys(gene_names)),
                "Organism": (item.get("organism") or {}).get("scientificName"),
                "Length": (item.get("sequence") or {}).get("length"),
            }
        )

    return pd.DataFrame(
        records,
        columns=["Accession", "Protein Name", "Gene", "Organism", "Length"],
    )


def collect_protein_info(query, output=None, size=10):
    if not query:
        raise ValueError("Protein/gene name is required.")

    df = search_uniprot(query, size=size)

    if df.empty:
        raise ValueError(f"No UniProt records found for '{query}'.")

    if output is None:
        output = f"{safe_filename(query)}_proteins.csv"

    df.to_csv(output, index=False)

    print(f"\nQuery: {query}")
    print(f"Records retrieved: {len(df)}")
    print(f"Saved: {output}\n")
    print(df.to_string(index=False))
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Retrieve protein information from UniProt."
    )
    parser.add_argument("protein", nargs="?", help="Protein/gene name, e.g. TP53")
    parser.add_argument(
        "-o", "--output", help="CSV output path (default: <query>_proteins.csv)"
    )
    parser.add_argument(
        "-n", "--size", type=int, default=10, help="Maximum number of records"
    )
    args = parser.parse_args()

    protein = args.protein
    if not protein:
        protein = input("Enter protein/gene name (e.g. TP53): ").strip()

    if not protein:
        sys.exit("Please enter a protein/gene name.")

    try:
        collect_protein_info(protein, args.output, max(1, min(args.size, 500)))
    except requests.RequestException as exc:
        sys.exit(f"Network/API error: {exc}")
    except Exception as exc:
        sys.exit(f"Error: {exc}")


if __name__ == "__main__":
    main()


# -------------------------------
# Combined Mini Project Menu
# -------------------------------

def run_mini_project():
    print("=" * 65)
    print("BIOMEDICAL DATA COLLECTOR MINI PROJECT")
    print("Name: Tuba Durdana")
    print("=" * 65)
    print("1. OPTION 1 — Disease Target Data Collector")
    print("2. OPTION 2 — Protein Information Collector")
    print("3. Exit")
    print("=" * 65)

    choice = input("Select an option (1/2/3): ").strip()

    if choice == "1":
        disease = input("Enter disease name (e.g. Breast Cancer): ").strip()
        if not disease:
            print("Disease name cannot be empty.")
            return
        collect_disease_targets(
            disease_name=disease,
            output="disease_targets.csv"
        )

    elif choice == "2":
        protein = input("Enter protein/gene name (e.g. TP53): ").strip()
        if not protein:
            print("Protein/gene name cannot be empty.")
            return
        collect_protein_info(
            query=protein,
            output=f"{protein.strip().replace(' ', '_')}_proteins.csv",
            size=10
        )

    elif choice == "3":
        print("Project closed.")
    else:
        print("Invalid option. Please choose 1, 2, or 3.")


if __name__ == "__main__":
    run_mini_project()


# ------------------------------------------------
# PROJECT AUTHOR
# ------------------------------------------------
# Tuba Durdana
