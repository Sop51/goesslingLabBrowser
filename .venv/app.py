from flask import Flask, render_template, request, redirect, url_for, flash, session, jsonify
from forms import GeneForm, GroupForm, UploadFileForm, LoginForm, AnnotationForm
import os
from dotenv import load_dotenv
import scanpy as sc
import io
import base64
import matplotlib.pyplot as plt
import pandas as pd
import re
from numpy.f2py.symbolic import as_string

# load environment variables from .env file
load_dotenv()

# initialize the app
app = Flask(__name__)
# configure the secret key
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY', 'default_secret_key')

# credentials from environment variables
USERNAME = os.getenv('USERNAME')
PASSWORD = os.getenv('PASSWORD')

# define and configure the upload folder destination
UPLOAD_FOLDER = 'uploads/'
app.config['UPLOAD_FILES_DEST'] = UPLOAD_FOLDER

# create upload folder if it doesn't exist
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

# define a function to plot the gene UMAP
def plot_gene_umap(gene, file):
    # load the adata object
    adata = sc.read(file)

    # Check if the gene exists in the AnnData object
    if gene not in adata.var_names:
        return None  # Return None if the gene is not found

    # create the plot
    plt.figure(figsize=(30, 30))  # set the size of the plot
    sc.pl.umap(adata, color=gene, show=False)  # show must be false to avoid rendering immediately
    # save the plot to a BytesIO object
    img = io.BytesIO()
    plt.savefig(img, format='png')
    plt.close()  # close the plot to free memory
    img.seek(0)  # seek to the start of the object
    # encode the image in base64
    plot_url = base64.b64encode(img.getvalue()).decode('utf-8')
    # return the plot URL for rendering
    return f"data:image/png;base64,{plot_url}"  # Data URL format

# define a function to plot the group umap
def plot_group_umap(group, file):
    # load the adata object
    adata = sc.read(file)
    # create the plot
    plt.figure(figsize=(30, 30))
    sc.pl.umap(
        adata,
        color=group,
        legend_loc="on data",
        frameon=False,
        legend_fontsize=10,
        legend_fontoutline=2,
        show=False  # Ensure not to show the plot
    )
    # save the plot to a BytesIO object
    img = io.BytesIO()
    plt.savefig(img, format='png', bbox_inches='tight')  # Add bbox_inches to ensure no cutting off
    plt.close()  # close the plot to free memory
    img.seek(0)  # move to the start of the object
    # encode the image in base64
    plot_url = base64.b64encode(img.getvalue()).decode('utf-8')
    # return the plot URL for rendering
    return f"data:image/png;base64,{plot_url}"

# define a function to generate the differential expression table based on the annotation
def generate_table(annotation, file):
    # load the adata object
    adata = sc.read_h5ad(file)
    # get the df with the appropriate information
    filtered_de_df = pd.DataFrame(adata.varm['filtered_de'], index=adata.var_names)
    # initialize the df to hold the dataframes
    df_dict = {}
    # loop through columns in the df
    for col in filtered_de_df.columns:
        match = re.match(r'(\d+):(.+)', col)  # regex to match digits followed by a colon
        if match:
            prefix = match.group(1)  # get the numeric prefix

            # add column to the corresponding df in the dictionary
            if prefix not in df_dict:
                df_dict[prefix] = pd.DataFrame()  # initialize as an empty df if the prefix is not already a key
            df_dict[prefix][match.group(2)] = filtered_de_df[col]  # add the column to the df

    # extract only the wanted columns
    final_df = df_dict[annotation][['log2Mean', 'log2FC', 'percentage', 'percentage_fold_change', 'auroc', 'mwu_pval']]

    # Reset the index to include 'gene' as a regular column
    final_df.reset_index(inplace=True)
    final_df.rename(columns={'index': 'gene'}, inplace=True)

    return(final_df)


# create the login route
@app.route('/login', methods=['GET', 'POST'])
def login():
    # set the error to none
    error = None
    # initialize the login form
    form = LoginForm()
    # check validation
    if form.validate_on_submit():
        # extract the usernames and see if they equal the env variables
        if form.username.data != USERNAME or form.password.data != PASSWORD:
            error = 'Invalid Credentials. Please try again.'
        else:
            session['logged_in'] = True
            return redirect(url_for('home'))
    return render_template('login.html', form=form, error=error)

# create the index route
@app.route('/', methods=['GET'])
def index():
    return redirect(url_for('login'))

# create the home route
@app.route('/home', methods=['GET', 'POST'])
def home():
    # Ensure the user is logged in
    if not session.get('logged_in'):
        return redirect(url_for('login'))

    # Create form instances
    fileform = UploadFileForm()
    annotationform = AnnotationForm()
    groupform = GroupForm()
    geneform = GeneForm()

    # list existing files in the uploads directory
    existing_files = os.listdir(app.config['UPLOAD_FILES_DEST'])

    # Handle file selection from existing files
    if request.method == 'POST':
        if fileform.validate_on_submit():
            # Get the file value from the form
            file = fileform.file.data
            filepath = os.path.join(app.config['UPLOAD_FILES_DEST'], file.filename)
            file.save(filepath)
            flash(f'"{file.filename}" uploaded successfully', 'success')

            # Reset annotations and selected group in the session after uploading
            session['filepath'] = filepath
            adata = sc.read(filepath)
            obs_columns = adata.obs.columns.tolist()
            annotcols = adata.obs['leiden'].unique().tolist()
            session['ann_cols'] = annotcols
            session['obs_columns'] = obs_columns
            return redirect(url_for('home'))  # Redirect to avoid resubmission
        if 'existing_file' in request.form:
            selected_file = request.form['existing_file']
            if selected_file:
                session['filepath'] = 'uploads/' + selected_file
                flash(f'File "{selected_file}" loaded successfully.', 'success')

                # Load the adata object and update session data as before
                adata = sc.read(session['filepath'])
                obs_columns = adata.obs.columns.tolist()
                annotcols = adata.obs['leiden'].unique().tolist()
                session['ann_cols'] = annotcols
                session['obs_columns'] = obs_columns
                return redirect(url_for('home'))

    # Populate group and annotation forms based on session-stored names
    if 'obs_columns' in session:
        groupform.group.choices = [(col, col) for col in session['obs_columns']]
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]

    return render_template('index.html', fileform=fileform, annotationform=annotationform, groupform=groupform,
                           geneform=geneform, existing_files=existing_files)

@app.route('/plot_gene', methods=['POST'])
def plot_gene():
    groupform = GroupForm()
    geneform = GeneForm()
    fileform = UploadFileForm()
    annotationform = AnnotationForm()


    if 'obs_columns' in session:
        groupform.group.choices = [(col, col) for col in session['obs_columns']]
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]

    if 'selected_group' in session:
        groupform.group.data = session['selected_group']
    if 'annotation' in session:
        annotationform.annotation.data = session['annotation']

    if geneform.validate_on_submit():
        gene_of_interest = geneform.gene.data
        session['gene_of_interest'] = gene_of_interest
        gene_plot = plot_gene_umap(gene_of_interest, session['filepath'])

        if gene_plot is None:
            # If gene is not found, flash an error message
            return jsonify({"error": f"Gene '{gene_of_interest}' not found in the dataset."})

        # Return JSON response
        return jsonify(gene_plot=gene_plot)

    return jsonify(error="Invalid form submission"), 400  # Return error if form is invalid


@app.route('/plot_group', methods=['POST'])
def plot_group():
    groupform = GroupForm()
    geneform = GeneForm()
    fileform = UploadFileForm()
    annotationform = AnnotationForm()

    if 'obs_columns' in session:
        groupform.group.choices = [(col, col) for col in session['obs_columns']]
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]

    gene_plot = None
    if 'gene_plot' in request.form:
        gene_plot = request.form['gene_plot']

    if 'gene_of_interest' in session:
        geneform.gene.data = session['gene_of_interest']
    if 'annotation' in session:
        annotationform.annotation.data = session['annotation']

    if groupform.validate_on_submit():
        selected_group = groupform.group.data
        session['selected_group'] = selected_group
        group_plot = plot_group_umap(selected_group, session['filepath'])

        # Return JSON response
        return jsonify(group_plot=group_plot)

    return jsonify(error="Invalid form submission"), 400  # Return error if form is invalid

@app.route('/make_table', methods=['POST'])
def make_table():
    print("make_table route called")
    annotationform = AnnotationForm()
    print("Form data:", request.form)
    if 'ann_cols' in session:
        annotationform.annotation.choices = [(annotation, annotation) for annotation in session['ann_cols']]
    if annotationform.validate_on_submit():
        annotation = annotationform.annotation.data
        print(f"Annotation submitted: {annotation}")
        session['annotation'] = annotation
        # Instead of redirecting, return a JSON response
        return jsonify({'status': 'success', 'message': 'Annotation saved'})
    else:
        print("Form validation failed:", annotationform.errors)
        return jsonify({'status': 'error', 'message': 'Form validation failed'}), 400


@app.route('/api/table_data', methods=['GET'])
def table_data():
    # get the selected annotation from the session
    annotation = session.get('annotation')
    filepath = session.get('filepath')

    print(f"Selected annotation: {annotation}")

    if annotation and filepath:
        # generate table data based on the annotation
        table_data = generate_table(annotation, filepath)
        # convert df to list of dictionaries for JSON response
        json_data = table_data.to_dict(orient='records')
        # save the current table to the session
        return jsonify(json_data)
    else:
        return jsonify([])  # return an empty JSON response if no data

# define the logout route
@app.route('/logout')
def logout():
    session.pop('logged_in', None)
    session.pop('filepath', None)
    session.pop('gene_of_interest', None)
    session.pop('selected_group', None)
    session.pop('selected_gene', None)
    session.pop('gene_plot', None)
    session.pop('group_plot', None)
    session.pop('obs_columns', None)
    session.pop('annotationtable', None)
    session.pop('ann_cols', None)
    return redirect(url_for('login'))

if __name__ == '__main__':
    app.run(debug=True)
