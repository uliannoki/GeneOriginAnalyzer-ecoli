import pandas as pd

# Пути к файлам
input_file = "results/synteny/candidate_neighborhoods.tsv"
output_file = "results/visualizations/gene_neighborhoods.csv"

# Твои 5 кандидатов (из логов хакатона)
candidates = [
    "ecoli_Sakai__NP_310248.1",
    "ecoli_Sakai__NP_303550.1",
    "ecoli_Sakai__NP_308259.1",
    "ecoli_REL606__ACT39985.1",
    "ecoli_Sakai__NP_310689.1"
]

# Читаем данные
df = pd.read_csv(input_file, sep="\t")
df = df[df["candidate_protein_id"].isin(candidates)]

rows = []

for cand_id, group in df.groupby("candidate_protein_id"):
    # Короткое имя вида (Sakai, REL606) и самого гена
    species = group["sample"].iloc[0].replace("ecoli_", "")
    gene_id = cand_id.split("__")[1]
    
    # 1. Сам кандидат
    cand_row = group[group["is_candidate"] == True].iloc[0]
    rows.append({
        "candidate": gene_id,
        "species": species,
        "gene": gene_id,
        "start": cand_row["candidate_start"],
        "end": cand_row["candidate_end"],
        "strand": cand_row["candidate_strand"],
        "role": "candidate"
    })
    
    # 2. Левый сосед (relative_order == -1)
    left = group[group["relative_order"] == -1]
    if not left.empty:
        l = left.iloc[0]
        rows.append({
            "candidate": gene_id,
            "species": species,
            "gene": l["neighbor_gene"] if pd.notna(l["neighbor_gene"]) else l["neighbor_protein_id"],
            "start": l["neighbor_start"],
            "end": l["neighbor_end"],
            "strand": l["neighbor_strand"],
            "role": "flanking"
        })
        
    # 3. Правый сосед (relative_order == 1)
    right = group[group["relative_order"] == 1]
    if not right.empty:
        r = right.iloc[0]
        rows.append({
            "candidate": gene_id,
            "species": species,
            "gene": r["neighbor_gene"] if pd.notna(r["neighbor_gene"]) else r["neighbor_protein_id"],
            "start": r["neighbor_start"],
            "end": r["neighbor_end"],
            "strand": r["neighbor_strand"],
            "role": "flanking"
        })

out_df = pd.DataFrame(rows)
out_df.to_csv(output_file, index=False)
print(f"Сгенерировано {len(out_df)} строк в {output_file}")
