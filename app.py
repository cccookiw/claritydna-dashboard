import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt

# Page setup
st.set_page_config(page_title="ClarityDNA", layout="wide")

# Load data
@st.cache_data
def load_data():
    df = pd.read_csv(
        "variant_summary_small.txt.gz",
        sep="\t",
        compression="gzip",
        usecols=["GeneSymbol", "ClinicalSignificance", "Type", "PhenotypeList", "ReviewStatus"],
        low_memory=False
    )
    return df

df = load_data()

# Clean data
df = df.dropna(subset=["GeneSymbol", "ClinicalSignificance"])
df = df.rename(columns={
    "GeneSymbol": "Gene",
    "ClinicalSignificance": "Risk",
    "Type": "Variant Type",
    "PhenotypeList": "Condition",
    "ReviewStatus": "Confidence"
})

# Simplify risk
def simplify_risk(x):
    x = str(x).lower()
    if "pathogenic" in x:
        return "High"
    elif "uncertain significance" in x:
        return "Moderate"
    elif "benign" in x:
        return "Low"
    else:
        return "Other"

df["Risk"] = df["Risk"].apply(simplify_risk)
df = df[df["Risk"].isin(["High", "Moderate", "Low"])]
df = df.head(20000)

# Title
st.title("🧬 ClarityDNA")
st.subheader("A Patient-Friendly Genomic Report Dashboard")
st.markdown("Turn complex genetic findings into clear, understandable results for patients.")

# Metrics
high = (df["Risk"] == "High").sum()
moderate = (df["Risk"] == "Moderate").sum()
low = (df["Risk"] == "Low").sum()

c1, c2, c3 = st.columns(3)
c1.metric("🔴 High Risk Findings", high)
c2.metric("🟡 Moderate Risk Findings", moderate)
c3.metric("🟢 Low Risk Findings", low)

if high > 0:
    st.warning("Some findings may deserve follow-up with a healthcare provider.")
else:
    st.success("No high-risk findings detected in the displayed data.")

# Sidebar filters
st.sidebar.header("Filter Results")

risk_filter = st.sidebar.selectbox("Risk Level", ["All", "High", "Moderate", "Low"])
gene_search = st.sidebar.text_input("Search by Gene")

# Categories
df["Category"] = df["Condition"].fillna("Unknown").apply(
    lambda x: "Cancer" if "cancer" in str(x).lower()
    else "Heart" if "cardio" in str(x).lower() or "heart" in str(x).lower()
    else "Neurological" if "alzheimer" in str(x).lower() or "neuro" in str(x).lower()
    else "Other"
)

category_filter = st.sidebar.selectbox(
    "Condition Category",
    ["All"] + sorted(df["Category"].dropna().unique().tolist())
)

# Apply filters
filtered_df = df.copy()

if risk_filter != "All":
    filtered_df = filtered_df[filtered_df["Risk"] == risk_filter]

if gene_search:
    filtered_df = filtered_df[
        filtered_df["Gene"].str.contains(gene_search, case=False, na=False)
    ]

if category_filter != "All":
    filtered_df = filtered_df[filtered_df["Category"] == category_filter]

# Layout
left, right = st.columns([2, 1])

# LEFT SIDE (Details)
with left:
    st.subheader("Variant Explorer")

    selected_gene = st.selectbox(
        "Select a gene to explore",
        ["Choose a gene"] + sorted(filtered_df["Gene"].dropna().unique().tolist())
    )

    if selected_gene != "Choose a gene":
        gene_data = filtered_df[filtered_df["Gene"] == selected_gene]

        for _, row in gene_data.head(5).iterrows():
            with st.container():
                st.markdown(f"### {row['Gene']}")
                st.write(f"**Risk Level:** {row['Risk']}")
                st.write(f"**Variant Type:** {row['Variant Type']}")
                st.write(f"**Condition:** {row['Condition']}")
                st.write(f"**Confidence:** {row['Confidence']}")

                if row["Risk"] == "High":
                    meaning = "This finding may be associated with a higher health risk."
                    action = "Consider discussing screening or follow-up with a doctor or genetic counselor."
                    worry = "This is not a diagnosis, but it may deserve attention."
                elif row["Risk"] == "Moderate":
                    meaning = "This finding has uncertain meaning and may need further interpretation."
                    action = "Ask a provider whether more review or testing is needed."
                    worry = "This does not mean disease is present."
                else:
                    meaning = "This finding is not usually expected to cause harm."
                    action = "No immediate action is typically needed."
                    worry = "This result is generally reassuring."

                st.write(f"📌 **What this means:** {meaning}")
                st.write(f"💡 **What to do next:** {action}")
                st.write(f"🧠 **Should you panic?** {worry}")
                st.markdown("---")
    else:
        st.info("Select a gene above to view patient-friendly details.")

# RIGHT SIDE (Charts)
with right:
    st.subheader("Risk Overview")

    risk_counts = filtered_df["Risk"].value_counts()

    fig, ax = plt.subplots()
    risk_counts.plot(kind="bar", ax=ax)
    ax.set_xlabel("Risk Level")
    ax.set_ylabel("Count")
    ax.set_title("Risk Distribution")
    st.pyplot(fig)

    st.subheader("Results Snapshot")
    st.write(f"Showing **{len(filtered_df):,}** filtered variants.")

# Info section
st.subheader("Understand Your Results")
st.markdown("""
- A **gene** is a section of DNA that helps control how your body works.
- A **variant** is a change in a gene.
- A **high-risk result does not mean you will definitely get a disease**.
- These results are meant to support conversations with healthcare professionals.
""")

# Disclaimer
st.subheader("Important Notice")
st.error("""
This dashboard is for education and communication purposes only.
It does not replace professional medical advice, diagnosis, or treatment.
Genetic findings should always be interpreted by a qualified healthcare professional.
""")
