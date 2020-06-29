
import textwrap
project = 'ZERO'
copyright = '2020, GD and SS'
author = 'GD and SS'

# The full version, including alpha/beta/rc tags
release = '1.0.0'
extensions = [
    # there may be others here already, e.g. 'sphinx.ext.mathjax'
    'breathe',
    'exhale',
    'recommonmark',
    'sphinx.ext.todo'
]



# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']
todo_include_todos = True
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_contex ={
    "display_github": True,
    "github_user": "ds4dm",
    "github_repo": "ZERO",
}
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background' : '#333131',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 5,
    'includehidden': True,
    'titles_only': False
}
html_logo = "support_files/zero_white.png"
# Setup the breathe extension
breathe_projects = {"My Project": "./doxyoutput/xml"}
breathe_default_project = "My Project"
breathe_default_members = ('members', 'undoc-members' , 'protected-members', 'private-members')


# Setup the exhale extension
exhale_args = {
    "containmentFolder":     "./api",
    "rootFileName":          "library_root.rst",
    "rootFileTitle":         "Library API",
    "doxygenStripFromPath":  "..",
    "createTreeView":        True,
    "exhaleExecutesDoxygen": True,
    "exhaleDoxygenStdin":    textwrap.dedent('''
        INPUT                   = ../README.md  ../src/ ../include
        IMAGE_PATH              = ./
        EXCLUDE_SYMBOLS         = arma::* rapidjson::*
        DOT_IMAGE_FORMAT        = png
        RECURSIVE               = YES
        EXTRACT_ALL             = YES
        EXTRACT_PRIVATE         = YES
        EXTERNAL_PAGES          = YES
        FILE_PATTERNS           = *.c *.h *.cpp
    ''')
}

# Tell sphinx what the primary language being documented is.
primary_domain = 'cpp'

# Tell sphinx what the pygments highlight language should be.
highlight_language = 'cpp'