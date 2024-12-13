# import the required packages
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileRequired, FileAllowed
from wtforms import StringField, SubmitField, SelectField, FileField, PasswordField, RadioField
from wtforms.validators import DataRequired, EqualTo, Length

# create a class for the login information
class LoginForm(FlaskForm):
    username = StringField('Username', validators=[DataRequired()])
    password = PasswordField('Password', validators=[DataRequired()])
    submit = SubmitField('Log In')

# create a class for the gene submission
class GeneForm(FlaskForm):
    gene = StringField('Gene of Interest', validators=[DataRequired()])
    submit = SubmitField('Submit')

# create a class for the group submission
class GroupForm(FlaskForm):
    group = SelectField('Group of Interest', choices=[], validators=[DataRequired()])
    submit = SubmitField('Submit')

# create a class for the file upload
class UploadFileForm(FlaskForm):
    file = FileField('Upload an h5Seurat file', validators=[
        FileRequired()
    ])
    submit = SubmitField('Upload')

# create a form to select an annotation for the table
class AnnotationForm(FlaskForm):
    annotation = SelectField('Annotation/Cluster of Interest', choices=[], validators=[DataRequired()])
    submit = SubmitField('Submit')

# create a form to select a gene for the subplot
class GeneSubplotForm(FlaskForm):
    gene = StringField('Gene of Interest', validators=[DataRequired()])
    submit = SubmitField('Submit')

# create a form to select the cell type of interest
class GroupSubplotForm(FlaskForm):
    group = SelectField('Cell Type of Interest', choices=[], validators=[DataRequired()])
    submit = SubmitField('Submit')

# create a form to display the timepoint plot
class TimepointForm(FlaskForm):
    timepoint = RadioField('View Timepoint Information', choices=[('Yes', 'yes'), ('No', 'no')], validators=[DataRequired()])
    submit = SubmitField('Submit')

class SubClusterForm(FlaskForm):
    cluster = SelectField('Cluster of Interest', choices=[], validators=[DataRequired()])
    submit = SubmitField('Submit')