from importlib.metadata import metadata
from flask import Flask, render_template, request, redirect, url_for, flash, session, jsonify
from forms import GeneForm, GroupForm, UploadFileForm, LoginForm, AnnotationForm
import os
from dotenv import load_dotenv
import io
import base64
import matplotlib.pyplot as plt
import pandas as pd
import re
from numpy.f2py.symbolic import as_string
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import tempfile

# load environment variables from .env file
load_dotenv()

# initialize the app
app = Flask(__name__)
# configure the secret key
app.config['SECRET_KEY'] = os.getenv('SECRET_KEY', 'default_secret_key')
app.config['WTF_CSRF_ENABLED'] = False

# credentials from environment variables
USERNAME = os.getenv('USERNAME')
PASSWORD = os.getenv('PASSWORD')

# define and configure the upload folder destination
UPLOAD_FOLDER = 'uploads/'
app.config['UPLOAD_FILES_DEST'] = UPLOAD_FOLDER

# create upload folder if it doesn't exist
if not os.path.exists(UPLOAD_FOLDER):
    os.makedirs(UPLOAD_FOLDER)

# create a global variable to hold the loaded seurat object
seurat_obj = None

# define a function to plot the gene UMAP
def plot_gene_umap(gene, file):
    pandas2ri.activate()
    # Load the required libraries in R
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    # load in the seurat obj
    global seurat_obj
    # Check if the gene exists in the Seurat object
    gene_check = ro.r(f'"{gene}" %in% rownames(seurat_obj)')
    if not gene_check[0]:
        return jsonify({"error": f"Gene '{gene}' not found in the dataset."})
    # Create the UMAP plot for the gene using Seurat's FeaturePlot
    ro.r(f'gene_plot <- FeaturePlot(seurat_obj, features = "{gene}", order=TRUE, min.cutoff = "q10", repel = TRUE)')
    # Create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # Use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # Print the plot to the device
        ro.r('print(gene_plot)')
        # Finalize the plot
        ro.r('dev.off()')
    # Now read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
    # Remove the temporary file after processing
    os.remove(temp_filename)
    # Return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"

# define a function to plot the group umap
def plot_group_umap(group, file):
    pandas2ri.activate()
    # Load the Seurat library in R
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')
    ro.r('library(ggplot2)')
    # load in the seurat obj
    global seurat_obj
    # Check if the group exists in the Seurat object metadata
    group_check = ro.r(f'"{group}" %in% colnames(seurat_obj@meta.data)')
    if not group_check[0]:
        return jsonify({"error": f"Group '{group}' not found in the dataset."})
    # Create the UMAP plot for the group using Seurat's DimPlot
    ro.r(
        f'group_plot <- DimPlot(seurat_obj, reduction = "umap", group.by = "{group}", label = TRUE, label.size = 5, repel = TRUE)')
    # Create a temporary file to save the plot
    with tempfile.NamedTemporaryFile(delete=False, suffix='.png') as temp_file:
        temp_filename = temp_file.name
        # Use the file name to create the plot
        ro.r(f'png("{temp_filename}", width = 800, height = 800, res = 150)')
        # Print the plot to the device
        ro.r('print(group_plot)')
        # Finalize the plot
        ro.r('dev.off()')
    # Now read the plot file and convert to base64
    with open(temp_filename, "rb") as img_file:
        img_data = img_file.read()
        plot_url = base64.b64encode(img_data).decode('utf-8')
    # Remove the temporary file after processing
    os.remove(temp_filename)
    # Return the plot as a base64-encoded string to embed in HTML
    return f"data:image/png;base64,{plot_url}"

# define a function to generate the differential expression table based on the annotation
def generate_table(annotation, file):
    pandas2ri.activate()
    # load the required seurat library
    ro.r('library(Seurat)')
    ro.r('library(SeuratDisk)')

    # load in the seurat obj
    global seurat_obj

    ro.r(f'''
    markers <- FindMarkers(seurat_obj, ident.1 = "{annotation}", group.by = "cell.type.9.long")
    markers <- as.data.frame(markers)
    markers$gene <- rownames(markers)  # Add row names as a new column
    rownames(markers) <- NULL
    markers <- markers[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
    ''')

    # Convert result into a pandas dataframe
    with (ro.default_converter + pandas2ri.converter).context():
        markers_df = ro.conversion.get_conversion().rpy2py(ro.r('markers'))
        markers_df.columns = markers_df.columns.str.strip()
        markers_df.columns = ['gene', 'avg_log2FC', 'pct1', 'pct2', 'p_val_adj']
        print("Markers DataFrame:", markers_df)
    return markers_df


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

    # define the global seurat obj
    global seurat_obj

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

            # load seurat
            ro.r('library(Seurat)')
            ro.r('library(SeuratDisk)')
            ro.r(f"seurat_obj <- LoadH5Seurat('{filepath}')")
            # store in the global variable
            seurat_obj = ro.r('seurat_obj')

            metadata_columns = []
            annotcols = []

            # get the col names for the umap plot
            with localconverter(pandas2ri.converter):
                meta_data_cols = ro.r(f'colnames(seurat_obj@meta.data)')
                metadata_columns = [str(col) for col in meta_data_cols]

            # get the annotation names for the table
            with localconverter(pandas2ri.converter):
                unique_values = ro.r(f'levels(seurat_obj@meta.data[["cell.type.9.long"]])')
                annotcols = [str(val) for val in unique_values]

            session['ann_cols'] = annotcols
            session['obs_columns'] = metadata_columns
            return redirect(url_for('home'))  # Redirect to avoid resubmission
        if 'existing_file' in request.form:
            selected_file = request.form['existing_file']
            if selected_file:
                session['filepath'] = 'uploads/' + selected_file
                filepath = session['filepath']
                flash(f'File "{selected_file}" loaded successfully.', 'success')

                # load seurat
                ro.r('library(Seurat)')
                ro.r('library(SeuratDisk)')
                ro.r(f"seurat_obj <- LoadH5Seurat('{filepath}')")
                # store in the global variable
                seurat_obj = ro.r('seurat_obj')

                print("Seurat object loaded:", seurat_obj)

                metadata_columns = []
                annotcols = []

                # get the col names for the umap plot
                with localconverter(pandas2ri.converter):
                    meta_data_cols = ro.r(f'colnames(seurat_obj@meta.data)')
                    metadata_columns = [str(col) for col in meta_data_cols]

                # get the annotation names for the table
                with localconverter(pandas2ri.converter):
                    unique_values = ro.r(f'levels(seurat_obj@meta.data[["cell.type.9.long"]])')
                    annotcols = [str(val) for val in unique_values]
                    print(unique_values)

                session['ann_cols'] = annotcols
                session['obs_columns'] = metadata_columns
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
    app.run(debug=True, threaded=False, use_reloader=False)
