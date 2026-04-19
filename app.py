import streamlit as st #creates interactive website
import pandas as pd #loads the dataset
import matplotlib.pyplot as plt #creates the visuals

# sets up the page, adds title, and a wide layout
st.set_page_config(page_title="ClarityDNA", layout="wide")

@st.cache_data #makes sure the data saves in memory so that the app is faster and we don't have to reload everytime
def load_data():
    df = pd.read_csv(
        "variant_summary_small.txt.gz", #genomic dataset we are using
        sep="\t", #separates file by tabs
        compression="gzip", #unizips compressed file
        usecols=["GeneSymbol", "ClinicalSignificance", "Type", "PhenotypeList", "ReviewStatus"], #keeps the important columns that we r using
        low_memory=False #we have a large data set so it helps not get errors
    )
    return df
#runs the function to load the data
df = load_data()

# Cleans data
df = df.dropna(subset=["GeneSymbol", "ClinicalSignificance"]) #removes rows with missing important info
df = df.rename(columns={
    "GeneSymbol": "Gene", #gene name
    "ClinicalSignificance": "Risk",  #risk level
    "Type": "Variant Type", #mutation type
    "PhenotypeList": "Condition", #condition
    "ReviewStatus": "Confidence" #reliability of result
})

# turns the complicated medical information into simple terms
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

df["Risk"] = df["Risk"].apply(simplify_risk) #applies the function to all the rows in the risk column
df = df[df["Risk"].isin(["High", "Moderate", "Low"])] #keeps main categories
df = df.head(20000) #limits size of the dataset

# Title and description of our interactive dashboard
st.title("🧬 ClarityDNA")
st.subheader("A Patient-Friendly Genomic Report Dashboard")
st.markdown("Turn complex genetic findings into clear, understandable results for patients.")

# counts the metrics in each risk category
high = (df["Risk"] == "High").sum()
moderate = (df["Risk"] == "Moderate").sum()
low = (df["Risk"] == "Low").sum()

c1, c2, c3 = st.columns(3) #creates 3 columns and shows the numbers side by side
#shows the counts with labels
c1.metric("🔴 High Risk Findings", high)
c2.metric("🟡 Moderate Risk Findings", moderate)
c3.metric("🟢 Low Risk Findings", low)
#creates warning for high risk results
if high > 0:
    st.warning("Some findings may deserve follow-up with a healthcare provider.")
else:
    st.success("No high-risk findings detected in the displayed data.")

# Sidebar filters
st.sidebar.header("Filter Results")
#dropdown to filter by the risk levels
risk_filter = st.sidebar.selectbox("Risk Level", ["All", "High", "Moderate", "Low"])
gene_search = st.sidebar.text_input("Search by Gene") #to search the specific gene name

# creates categories to be easier for users
df["Category"] = df["Condition"].fillna("Unknown").apply(
    lambda x: "Cancer" if "cancer" in str(x).lower()
    else "Heart" if "cardio" in str(x).lower() or "heart" in str(x).lower()
    else "Neurological" if "alzheimer" in str(x).lower() or "neuro" in str(x).lower()
    else "Other"
)
#dropdown to filter categories
category_filter = st.sidebar.selectbox(
    "Condition Category",
    ["All"] + sorted(df["Category"].dropna().unique().tolist())
)

# Applies filters
filtered_df = df.copy() #makes copy to not mess up original data

if risk_filter != "All":    #applies risk filter
    filtered_df = filtered_df[filtered_df["Risk"] == risk_filter]

if gene_search:     #applies the gene search filter
    filtered_df = filtered_df[
        filtered_df["Gene"].str.contains(gene_search, case=False, na=False)
    ]

if category_filter != "All":   #applies category filter
    filtered_df = filtered_df[filtered_df["Category"] == category_filter]

# Layout to make left side bigger and right side smaller
left, right = st.columns([2, 1])

# left side which are the details
with left:
    st.subheader("Variant Explorer")

    selected_gene = st.selectbox(  #dropdown to choose a gene
        "Select a gene to explore",
        ["Choose a gene"] + sorted(filtered_df["Gene"].dropna().unique().tolist())
    )

    if selected_gene != "Choose a gene": #shows the results only if the gene is selected
        gene_data = filtered_df[filtered_df["Gene"] == selected_gene]

        for _, row in gene_data.head(5).iterrows():  #loops through first 5 results
            with st.container():
                #displays information cleanly
                st.markdown(f"### {row['Gene']}")
                st.write(f"**Risk Level:** {row['Risk']}")
                st.write(f"**Variant Type:** {row['Variant Type']}")
                st.write(f"**Condition:** {row['Condition']}")
                st.write(f"**Confidence:** {row['Confidence']}")
                #translates the complicated medical information into simple terms
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
                #gives explanation to user
                st.write(f"📌 **What this means:** {meaning}")
                st.write(f"💡 **What to do next:** {action}")
                st.write(f"🧠 **Should you panic?** {worry}")
                st.markdown("---")
    else:
        st.info("Select a gene above to view patient-friendly details.")

# right side which are the charts
with right:
    st.subheader("Risk Overview")
    #counts how many of each risk types
    risk_counts = filtered_df["Risk"].value_counts()
    #creates bar chart
    fig, ax = plt.subplots()
    risk_counts.plot(kind="bar", ax=ax)
    #labels the chart
    ax.set_xlabel("Risk Level")
    ax.set_ylabel("Count")
    ax.set_title("Risk Distribution")
    st.pyplot(fig) #displays the chart in streamlit

    st.subheader("Results Snapshot")
    st.write(f"Showing **{len(filtered_df):,}** filtered variants.")

# information section
st.subheader("Understand Your Results")
st.markdown("""
- A **gene** is a section of DNA that helps control how your body works.
- A **variant** is a change in a gene.
- A **high-risk result does not mean you will definitely get a disease**.
- These results are meant to support conversations with healthcare professionals.
""")
