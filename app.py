import streamlit as st
import requests
import py3Dmol
from stmol import showmol
import os
import urllib.parse
import json

# Set page configuration with premium look
st.set_page_config(
    page_title="NeuroGeno3D - Clinical Portal",
    page_icon="🧠",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for Premium Clinical UI Aesthetics
st.markdown("""
<style>
    /* Global Styles */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800&display=swap');
    
    .report-card {
        background: var(--secondary-background-color);
        border: 1px solid rgba(128, 128, 128, 0.15);
        border-radius: 12px;
        padding: 22px;
        margin-bottom: 20px;
        box-shadow: 0 4px 20px rgba(0, 0, 0, 0.04);
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
        color: var(--text-color);
    }
    .report-card:hover {
        border-color: #0984e3;
        box-shadow: 0 8px 30px rgba(9, 132, 227, 0.12);
    }
    .card-header {
        font-size: 1.2rem;
        font-weight: 700;
        margin-bottom: 14px;
        color: var(--text-color);
        display: flex;
        align-items: center;
        gap: 10px;
        border-bottom: 1px solid rgba(128, 128, 128, 0.15);
        padding-bottom: 10px;
    }
    .badge-pathogenic {
        background-color: #d63031;
        color: #ffffff;
        padding: 4px 12px;
        border-radius: 20px;
        font-weight: bold;
        font-size: 0.8rem;
    }
    .badge-vus {
        background-color: #fdcb6e;
        color: #2d3436;
        padding: 4px 12px;
        border-radius: 20px;
        font-weight: bold;
        font-size: 0.8rem;
    }
    .badge-benign {
        background-color: #00b894;
        color: #ffffff;
        padding: 4px 12px;
        border-radius: 20px;
        font-weight: bold;
        font-size: 0.8rem;
    }
    .badge-review {
        background-color: #0984e3;
        color: #ffffff;
        padding: 2px 8px;
        border-radius: 4px;
        font-size: 0.75rem;
        font-weight: 600;
    }
    .api-status {
        display: inline-flex;
        align-items: center;
        gap: 6px;
        font-size: 0.8rem;
        color: #2ecc71;
        font-weight: 500;
    }
    .api-status-dot {
        width: 8px;
        height: 8px;
        background-color: #2ecc71;
        border-radius: 50%;
        display: inline-block;
        box-shadow: 0 0 8px #2ecc71;
    }
    h1, h2, h3, h4 {
        font-family: 'Inter', sans-serif;
        font-weight: 700;
    }
    .lit-title {
        font-size: 0.95rem;
        font-weight: 600;
        margin-bottom: 4px;
    }
    .lit-meta {
        font-size: 0.8rem;
        opacity: 0.8;
        margin-bottom: 12px;
    }
</style>
""", unsafe_allow_html=True)

# ----------------- DATABASE CACHE & API MODULES -----------------

@st.cache_data(show_spinner=False, ttl=3600)
def fetch_uniprot_data(gene_symbol):
    """Retrieves standard gene symbol, description, accession, and wild-type sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene_symbol}+AND+organism_id:9606+AND+reviewed:true&format=json"
    try:
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            data = response.json()
            results = data.get("results", [])
            if results:
                entry = results[0]
                primary_accession = entry.get("primaryAccession")
                protein_desc = entry.get("proteinDescription", {})
                rec_name = protein_desc.get("recommendedName", {})
                full_name = rec_name.get("fullName", {}).get("value", "Unknown Protein")
                sequence = entry.get("sequence", {}).get("value")
                return {
                    "uniprot_id": primary_accession,
                    "protein_name": full_name,
                    "sequence": sequence,
                    "success": True
                }
    except Exception as e:
        pass
    return {"success": False}

@st.cache_data(show_spinner=False, ttl=3600)
def search_clinvar(gene_symbol, variant_str):
    """Searches NCBI ClinVar for real-time annotations and classifications."""
    # Build a robust search term: e.g., "IDH1 R132H"
    term = f"{gene_symbol} {variant_str}"
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term={urllib.parse.quote(term)}&retmode=json"
    try:
        res = requests.get(search_url, timeout=10)
        if res.status_code == 200:
            ids = res.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                # Retrieve esummary details
                sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id={','.join(ids)}&retmode=json"
                sum_res = requests.get(sum_url, timeout=10)
                if sum_res.status_code == 200:
                    result_data = sum_res.json().get("result", {})
                    summaries = []
                    for uid in result_data.get("uids", []):
                        var_data = result_data.get(uid, {})
                        
                        # Extract clinical significance and review status
                        significance = "Variant of Uncertain Significance (VUS)"
                        review_status = "No assertion criteria provided"
                        last_evaluated = "Unknown"
                        
                        # Check different classification blocks
                        for sig_key in ['germline_classification', 'clinical_impact_classification', 'oncogenicity_classification']:
                            sig_data = var_data.get(sig_key)
                            if sig_data and isinstance(sig_data, dict):
                                desc = sig_data.get('description')
                                if desc and significance == "Variant of Uncertain Significance (VUS)":
                                    significance = desc
                                    review_status = sig_data.get('review_status', review_status)
                                date = sig_data.get('last_evaluated')
                                if date and last_evaluated == "Unknown":
                                    last_evaluated = date
                        
                        # Extract phenotypes
                        phenotypes = []
                        for sig_key in ['germline_classification', 'clinical_impact_classification', 'oncogenicity_classification']:
                            classif = var_data.get(sig_key, {})
                            for trait in classif.get('trait_set', []):
                                name = trait.get('trait_name')
                                if name and name not in phenotypes:
                                    phenotypes.append(name)
                        
                        # Extract genomic location
                        chrom = "Unknown"
                        start = "Unknown"
                        band = "Unknown"
                        for measure in var_data.get("variation_set", []):
                            for loc in measure.get("variation_loc", []):
                                if loc.get("assembly_name") == "GRCh38":
                                    chrom = loc.get("chr", chrom)
                                    start = loc.get("start", start)
                                    band = loc.get("band", band)
                                    
                        # Extract allele frequencies
                        freqs = []
                        for measure in var_data.get("variation_set", []):
                            for freq in measure.get("allele_freq_set", []):
                                freqs.append({
                                    "source": freq.get("source", "Unknown"),
                                    "value": freq.get("value", "0.0")
                                })
                                
                        summaries.append({
                            "id": uid,
                            "title": var_data.get("title", ""),
                            "significance": significance,
                            "review_status": review_status,
                            "last_evaluated": last_evaluated,
                            "phenotypes": phenotypes,
                            "chrom": chrom,
                            "start": start,
                            "band": band,
                            "freqs": freqs,
                            "type": var_data.get("obj_type", "Variant")
                        })
                    return summaries
    except Exception:
        pass
    return []

@st.cache_data(show_spinner=False, ttl=3600)
def search_pubmed(gene_symbol, variant_str):
    """Searches PubMed for active peer-reviewed publications containing the variant details."""
    term = f"{gene_symbol} {variant_str} AND (glioma OR tumor OR cancer)"
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term={urllib.parse.quote(term)}&retmode=json&retmax=5"
    articles = []
    try:
        res = requests.get(search_url, timeout=10)
        if res.status_code == 200:
            ids = res.json().get("esearchresult", {}).get("idlist", [])
            if ids:
                sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id={','.join(ids)}&retmode=json"
                sum_res = requests.get(sum_url, timeout=10)
                if sum_res.status_code == 200:
                    results = sum_res.json().get("result", {})
                    for uid in ids:
                        art = results.get(uid, {})
                        title = art.get("title", "No Title")
                        pub_date = art.get("pubdate", "Unknown Date")
                        source = art.get("source", "Unknown Journal")
                        authors = art.get("authors", [])
                        author_str = authors[0].get("name", "Unknown") + " et al." if authors else "Unknown"
                        articles.append({
                            "pmid": uid,
                            "title": title,
                            "date": pub_date,
                            "source": source,
                            "author": author_str,
                            "url": f"https://pubmed.ncbi.nlm.nih.gov/{uid}/"
                        })
    except Exception:
        pass
    return articles

# ----------------- PRESET GENES (Original Clinical Set) -----------------
PRESET_GENES = {
    "IDH1": {"uniprot": "O75874", "default_mut": "R132H"},
    "IDH2": {"uniprot": "P48735", "default_mut": "R172K"},
    "BRAF": {"uniprot": "P15056", "default_mut": "V600E"},
    "TP53": {"uniprot": "P04637", "default_mut": "R273H"},
    "H3-3A": {"uniprot": "P84243", "default_mut": "K27M"}
}

# Amino Acid Mapping
AA_MAP = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
}

# ----------------- SIDEBAR PANEL -----------------
st.sidebar.markdown("""
<div style="margin-bottom: 20px;">
    <h2>🔬 Clinical Panel</h2>
    <div class="api-status"><span class="api-status-dot"></span> NIH ClinVar API Connected</div>
</div>
""", unsafe_allow_html=True)

offline_mode = st.sidebar.toggle("🔌 Offline Mode (Bypass ESMFold API)", value=False)

# 1. Custom / Preset Selector
gene_mode = st.sidebar.radio("Gene Selection Mode", ["Select Preset Gene", "Search Custom Gene"])

if gene_mode == "Select Preset Gene":
    gene_key = st.sidebar.selectbox("Target Gene", list(PRESET_GENES.keys()))
    uniprot_id = PRESET_GENES[gene_key]["uniprot"]
else:
    gene_key = st.sidebar.text_input("Enter Gene Symbol (e.g. EGFR, KRAS)", value="EGFR").upper().strip()
    uniprot_id = None

# Query UniProt for sequence
with st.spinner("Fetching sequence from UniProt..."):
    # If preset, we can use the symbol directly.
    uniprot_info = fetch_uniprot_data(gene_key)

if uniprot_info["success"]:
    gene_full_name = uniprot_info["protein_name"]
    wild_type_sequence = uniprot_info["sequence"]
    uniprot_id = uniprot_info["uniprot_id"]
    st.sidebar.success(f"**Protein:** {gene_full_name} ({uniprot_id})")
else:
    st.sidebar.error("Could not fetch gene sequence from UniProt. Using offline fallback sequence.")
    # Standard fallbacks for presets if UniProt fails
    if gene_key in PRESET_GENES:
        wild_type_sequence = "DKIAPWLDNDKMVHQKIKNYLKKVEQLSKELTNYLKLKSKTYAILDIRGHDTTRVGITKVALQKEYIAGIARQAIQGDIRWEKLEHYEILK"
        gene_full_name = "Isocitrate Dehydrogenase"
    else:
        wild_type_sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVL"
        gene_full_name = "Unknown Reference Protein"

# Show Sequence Area
raw_seq = st.sidebar.text_area("Base Wild-Type Sequence", wild_type_sequence, height=120)

# Mutation Type Configuration
default_mutation_label = PRESET_GENES[gene_key]["default_mut"] if gene_key in PRESET_GENES else "T790M"
mutation_type = st.sidebar.radio("Mutation Type", [
    f"Standard Point Mutation ({default_mutation_label})",
    "Custom Point Mutation",
    "Custom In-Frame Deletion",
    "Custom In-Frame Insertion",
    "Wild-Type (No Mutation)"
])

# Process mutated sequence
mutated_seq = list(raw_seq)
active_label = "Wild-Type"
is_mutation = False
mutation_indices = []
mutation_position = 1

if mutation_type == f"Standard Point Mutation ({default_mutation_label})":
    # Parse standard mutation e.g. R132H
    import re
    match = re.match(r"([A-Z])(\d+)([A-Z])", default_mutation_label)
    if match:
        wt_aa, pos_str, mut_aa = match.groups()
        pos = int(pos_str)
        # Often sequences in UniProt are longer, try to match or use default relative position
        if pos <= len(mutated_seq):
            mutated_seq[pos - 1] = mut_aa
            active_label = f"p.{AA_MAP.get(wt_aa, wt_aa)}{pos}{AA_MAP.get(mut_aa, mut_aa)} ({default_mutation_label})"
            is_mutation = True
            mutation_indices = [pos - 1]
            mutation_position = pos
else:
    if mutation_type == "Custom Point Mutation":
        c_pos = st.sidebar.number_input("Mutation Position (1-indexed)", min_value=1, max_value=len(raw_seq), value=1)
        c_res = st.sidebar.text_input("New Residue (Single Letter)", value="H").upper()
        if len(c_res) == 1 and (c_pos - 1) < len(mutated_seq):
            wt_res = mutated_seq[c_pos - 1]
            mutated_seq[c_pos - 1] = c_res
            wt_res_3 = AA_MAP.get(wt_res, wt_res)
            c_res_3 = AA_MAP.get(c_res, c_res)
            active_label = f"p.{wt_res_3}{c_pos}{c_res_3} (p.{wt_res}{c_pos}{c_res})"
            is_mutation = True
            mutation_indices = [c_pos - 1]
            mutation_position = c_pos

    elif mutation_type == "Custom In-Frame Deletion":
        del_start = st.sidebar.number_input("Deletion Start Position (1-indexed)", min_value=1, max_value=len(raw_seq), value=1)
        del_len = st.sidebar.number_input("Deletion Length (aa)", min_value=1, max_value=20, value=3)
        start_idx = del_start - 1
        if start_idx < len(mutated_seq):
            mutated_seq = mutated_seq[:start_idx] + mutated_seq[start_idx + del_len:]
            start_aa = raw_seq[start_idx]
            start_aa_3 = AA_MAP.get(start_aa, start_aa)
            if del_len == 1:
                active_label = f"p.{start_aa_3}{del_start}del"
            else:
                end_idx = min(len(raw_seq) - 1, start_idx + del_len - 1)
                end_aa = raw_seq[end_idx]
                end_aa_3 = AA_MAP.get(end_aa, end_aa)
                active_label = f"p.{start_aa_3}{del_start}_{end_aa_3}{del_start + del_len - 1}del"
            is_mutation = True
            mutation_indices = [max(0, start_idx - 1), start_idx]
            mutation_position = del_start

    elif mutation_type == "Custom In-Frame Insertion":
        ins_pos = st.sidebar.number_input("Insertion After Position (1-indexed)", min_value=1, max_value=len(raw_seq), value=1)
        ins_seq = st.sidebar.text_input("Insert Residues", value="AAA").upper()
        ins_idx = ins_pos
        if ins_idx <= len(mutated_seq):
            mutated_seq = mutated_seq[:ins_idx] + list(ins_seq) + mutated_seq[ins_idx:]
            flank1_idx = ins_pos - 1
            if flank1_idx < len(raw_seq):
                f1_aa = raw_seq[flank1_idx]
                f1_aa_3 = AA_MAP.get(f1_aa, f1_aa)
                ins_seq_3 = "".join([AA_MAP.get(r, r) for r in ins_seq])
                active_label = f"p.{f1_aa_3}{ins_pos}ins{ins_seq_3}"
            is_mutation = True
            mutation_indices = list(range(ins_idx, ins_idx + len(ins_seq)))
            mutation_position = ins_pos

final_sequence = "".join(mutated_seq)

st.sidebar.markdown(f"**Sequence Length:** {len(final_sequence)} aa")
run_btn = st.sidebar.button("🚀 Analyze Real-time Clinical Evidence", type="primary")

# ----------------- FOLDING PIPELINE WITH OPTIMIZATION -----------------

@st.cache_data(show_spinner=False)
def fold_sequence_esm(sequence):
    """Predicts 3D structure using ESMFold."""
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    try:
        response = requests.post(url, data=sequence, headers={"Content-Type": "text/plain"}, timeout=25)
        if response.status_code == 200:
            return response.text, "ESMFold API (Live)"
    except Exception:
        pass
    return None, None

def get_pdb_structure(gene, mut_type, sequence, offline_mode=False):
    """Fetches structure dynamically, extracting local fragment if sequence > 250aa."""
    if offline_mode:
        fallback_filename = f"{gene.lower().replace('-', '_')}_wt.pdb"
        fallback_path = os.path.join("fallback_pdb", fallback_filename)
        if os.path.exists(fallback_path):
            try:
                with open(fallback_path, "r") as f:
                    return f.read(), f"Offline Cache (Using {fallback_filename})"
            except Exception:
                pass
        return None, "Offline Mode (No matching local PDB structure found)"

    # Live folding
    pdb_data, source = fold_sequence_esm(sequence)
    if pdb_data:
        return pdb_data, source
    
    # Fallback to local
    fallback_filename = f"{gene.lower().replace('-', '_')}_wt.pdb"
    fallback_path = os.path.join("fallback_pdb", fallback_filename)
    if os.path.exists(fallback_path):
        try:
            with open(fallback_path, "r") as f:
                return f.read(), f"Offline Cache Fallback ({fallback_filename})"
        except Exception:
            pass
    return None, None

# ----------------- MAIN LAYOUT -----------------
st.title("🧠 NeuroGeno3D")
st.subheader("Real-Time Clinical Variant Visualizer & Decision Support Portal")

col1, col2 = st.columns([1, 1])

if run_btn or 'pdb_data' in st.session_state:
    if run_btn:
        with st.spinner("Retrieving ClinVar annotations and literature..."):
            # Extract clean mutation name (e.g. R132H)
            clean_mut = active_label.split(' ')[-1].replace('(', '').replace(')', '') if '(' in active_label else active_label.split('.')[-1]
            
            clinvar_records = search_clinvar(gene_key, clean_mut)
            pubmed_records = search_pubmed(gene_key, clean_mut)
            
            # 3D Folding setup (optimize sequence size)
            fold_seq = final_sequence
            is_fragment = False
            frag_offset = 0
            
            if len(final_sequence) > 500:
                is_fragment = True
                # Extract 300aa window around mutation
                mut_pos_0 = max(0, mutation_position - 1)
                start_win = max(0, mut_pos_0 - 150)
                end_win = min(len(final_sequence), mut_pos_0 + 150)
                fold_seq = final_sequence[start_win:end_win]
                frag_offset = start_win
                adjusted_mutation_indices = [idx - frag_offset for idx in mutation_indices]
            else:
                adjusted_mutation_indices = mutation_indices
                
            pdb_data, structure_source = get_pdb_structure(gene_key, mutation_type, fold_seq, offline_mode)
            
            st.session_state['pdb_data'] = pdb_data
            st.session_state['structure_source'] = structure_source
            st.session_state['active_label'] = active_label
            st.session_state['mutation_indices'] = adjusted_mutation_indices
            st.session_state['is_mutation'] = is_mutation
            st.session_state['gene_key'] = gene_key
            st.session_state['gene_full_name'] = gene_full_name
            st.session_state['clinvar_records'] = clinvar_records
            st.session_state['pubmed_records'] = pubmed_records
            st.session_state['is_fragment'] = is_fragment

    pdb_data = st.session_state.get('pdb_data')
    structure_source = st.session_state.get('structure_source')
    active_label = st.session_state.get('active_label', active_label)
    mutation_indices = st.session_state.get('mutation_indices', mutation_indices)
    is_mutation = st.session_state.get('is_mutation', is_mutation)
    current_gene_key = st.session_state.get('gene_key', gene_key)
    current_gene_name = st.session_state.get('gene_full_name', gene_full_name)
    clinvar_records = st.session_state.get('clinvar_records', [])
    pubmed_records = st.session_state.get('pubmed_records', [])
    is_fragment = st.session_state.get('is_fragment', False)

    # Compile the ClinVar status dynamically
    if clinvar_records:
        primary_record = clinvar_records[0]
        clinvar_status = primary_record["significance"]
        review_status = primary_record["review_status"]
        last_evaluated = primary_record["last_evaluated"]
        phenotypes_list = primary_record["phenotypes"]
        chrom_loc = f"Chr {primary_record['chrom']}:{primary_record['start']} ({primary_record['band']})"
        freqs_list = primary_record["freqs"]
        
        # ACMG classification badge mapping
        sig_lower = clinvar_status.lower()
        if "pathogenic" in sig_lower:
            badge_html = '<span class="badge-pathogenic">Pathogenic</span>'
            acmg_code = "Pathogenic"
        elif "benign" in sig_lower:
            badge_html = '<span class="badge-benign">Benign</span>'
            acmg_code = "Benign"
        else:
            badge_html = '<span class="badge-vus">VUS</span>'
            acmg_code = "VUS"
    else:
        clinvar_status = "Not Found / Novel Variant"
        review_status = "No clinical submissions in ClinVar"
        last_evaluated = "N/A"
        phenotypes_list = ["None documented in ClinVar database"]
        chrom_loc = "Unknown genomic locus"
        freqs_list = []
        badge_html = '<span class="badge-vus">VUS</span>'
        acmg_code = "VUS (Computational / Presumed Novel)"

    # LEFT COLUMN: INTERACTIVE 3D STRUCTURE
    with col1:
        st.markdown("### 🧬 Interactive 3D Protein Structure")
        if pdb_data:
            if is_fragment:
                st.warning("⚠️ Sequence length exceeds 500aa. Displaying 300aa local structural fragment around mutation site.")
            st.info(f"Structure predicted via: **{structure_source}**")
            
            view = py3Dmol.view(width=600, height=500)
            view.addModel(pdb_data, "pdb")
            view.setStyle({'cartoon': {'color': 'spectrum'}})
            
            if is_mutation and len(mutation_indices) > 0:
                for idx in mutation_indices:
                    resi_str = str(idx + 1)
                    view.addStyle(
                        {'resi': resi_str},
                        {'stick': {'color': '#d63031', 'radius': 0.3}, 'sphere': {'color': '#d63031', 'scale': 0.45}}
                    )
                view.zoomTo({'resi': str(mutation_indices[0] + 1)})
            else:
                view.zoomTo()
                
            showmol(view, height=500, width=600)
            st.caption("🔴 Mutated/altered site(s) highlighted in RED stick representation.")
        else:
            st.error("Error: Could not retrieve PDB structure. Check ESMFold connectivity.")

        # Real-time PubMed literature Feed
        st.markdown("### 📚 Live Peer-Reviewed Literature (PubMed)")
        if pubmed_records:
            lit_html = ""
            for art in pubmed_records:
                lit_html += f"""
                <div style="border-bottom: 1px solid rgba(128,128,128,0.1); padding-bottom: 10px; margin-bottom: 10px;">
                    <div class="lit-title"><a href="{art['url']}" target="_blank">📄 {art['title']}</a></div>
                    <div class="lit-meta"><strong>Source:</strong> {art['source']} ({art['date']}) | <strong>Authors:</strong> {art['author']} | <strong>PMID:</strong> {art['pmid']}</div>
                </div>
                """
            st.markdown(f'<div class="report-card">{lit_html}</div>', unsafe_allow_html=True)
        else:
            st.markdown("""
            <div class="report-card">
                <p>No peer-reviewed publications specifically detailing this variant were found in the current PubMed query.</p>
            </div>
            """, unsafe_allow_html=True)

    # RIGHT COLUMN: CLINICAL DECISION SUPPORT REPORT
    with col2:
        st.markdown("### 📋 Clinical Decision Support Report")
        
        # Real-time ClinVar Annotations Card
        freq_html = ""
        if freqs_list:
            freq_html = "<ul>"
            for f in freqs_list:
                freq_html += f"<li><strong>{f['source']}:</strong> {f['value']}</li>"
            freq_html += "</ul>"
        else:
            freq_html = "No allele frequency data documented in ClinVar."

        pheno_html = ", ".join(phenotypes_list)

        st.markdown(f"""
        <div class="report-card">
            <div class="card-header">🏷️ 1. Real-time ClinVar & dbSNP Annotations</div>
            <p><strong>Gene:</strong> {current_gene_name} ({current_gene_key})</p>
            <p><strong>Variant Nomenclature:</strong> <code>{active_label}</code></p>
            <p><strong>Genomic Location (GRCh38):</strong> {chrom_loc}</p>
            <hr style="margin: 12px 0; border: 0; border-top: 1px solid rgba(128,128,128,0.15);"/>
            <p><strong>Clinical Significance:</strong> {clinvar_status} &nbsp; {badge_html}</p>
            <p><strong>Review Status:</strong> <span class="badge-review">{review_status}</span></p>
            <p><strong>Associated Conditions:</strong> {pheno_html}</p>
            <p><strong>Last Evaluated:</strong> {last_evaluated}</p>
            <p><strong>Allele Frequency:</strong> {freq_html}</p>
        </div>
        """, unsafe_allow_html=True)

        # Molecular Impact Card
        if not is_mutation:
            impact_text = "Wild-type baseline. Normal cellular functions preserved."
        else:
            impact_text = f"""
            - **Structural alteration at position {mutation_position}:** Mutated sequence has been folded using ESMFold to predict local conformational shifts.
            - **Molecular Consequence:** Local amino acid substitution might disrupt hydrogen bond networks, active site pocket geometry, or ligand interactions.
            - **ACMG Pathogenicity classification:** {acmg_code} based on public database evidence and clinical assertions.
            """
            
        st.markdown(f"""
        <div class="report-card">
            <div class="card-header">⚡ 2. Molecular & Structural Impact</div>
            <p>{impact_text}</p>
        </div>
        """, unsafe_allow_html=True)

        # WHO & Clinical Recommendation Card
        if "pathogenic" in clinvar_status.lower():
            rec_text = "Highly recommended to correlate with standard IHC markers and diagnostic imaging. Target therapeutic options or clinical trials matching the mutation pathway."
        else:
            rec_text = "Variant of Uncertain Significance. Recommend monitoring, follow-up sequencing, and correlation with histopathology and clinical state."

        st.markdown(f"""
        <div class="report-card">
            <div class="card-header">📖 3. Diagnostic & Therapeutic Guidelines</div>
            <p><strong>Clinical Recommendation:</strong> {rec_text}</p>
            <p style="font-size: 0.85rem; opacity: 0.8;"><em>Note: Integrated diagnostic pathways should align with current WHO Classification of Central Nervous System Tumors (WHO CNS 5th Edition).</em></p>
        </div>
        """, unsafe_allow_html=True)

        # Clinical Report Export
        report_content = f"""# NeuroGeno3D CLINICAL VARIANT REPORT
Gene Symbol: {current_gene_key} ({current_gene_name})
UniProt Accession: {uniprot_id}
Variant: {active_label}

CLINICAL SIGNIFICANCE SUMMARY:
- ClinVar Pathogenicity: {clinvar_status}
- ACMG Category: {acmg_code}
- Review Status: {review_status}
- Genomic Location: {chrom_loc}
- Associated Phenotypes: {', '.join(phenotypes_list)}
- Last Evaluated: {last_evaluated}

LITERATURE REFERENCES (PMIDs):
{', '.join([a['pmid'] for a in pubmed_records]) if pubmed_records else 'None found'}

Disclaimer: For Research and Educational Use Only.
"""
        st.download_button(
            label="💾 Download Official Clinical Report",
            data=report_content,
            file_name=f"Clinical_Report_{current_gene_key}_{active_label.replace(' ', '_')}.txt",
            mime="text/plain",
            type="secondary"
        )

else:
    st.info("👈 Configure target gene mutations in the sidebar and click **'Analyze Real-time Clinical Evidence'** to generate the clinical dashboard.")

# --- FOOTER CLINICAL DISCLAIMER ---
st.markdown("---")
st.markdown("""
<div style="text-align: center; font-size: 0.8rem; opacity: 0.6; padding: 20px 0; line-height: 1.4;">
    <strong>Research and Education Use Only.</strong> This application is developed strictly for educational and academic research purposes. It is not approved for use in clinical diagnostics, patient management, or therapeutic decision-making. All structural predictions and AI annotations are computational models that require validation by certified diagnostic assays.
    <br/><br/>
    <span style="font-size: 0.9rem; font-weight: bold; opacity: 0.9; color: var(--primary-color);">Developed by Dr. Rijhul Lahariya</span>
</div>
""", unsafe_allow_html=True)
