# A beginner’s guide to coding style

Created: November 13, 2024 3:24 PM

# Why style matters

Having a good coding style means that it makes it easier for others to understand, use, and modify your code. Here “others” including your future self - so even when you don’t expect sharing your code with anyone, taking some extra seconds (or minutes) to maintain good coding style and reflect on best practices could save you hours/days/months in the future, especially when you are working on a large and/or long-term project.

The purpose of this note is to give a short summary of some considerations and best practices, particularly in the context of astrophysics data analysis using python and jupyter notebooks and for researchers who have not received extensive formal training in software development.

For more extensive guides in various contexts, see:

- [google’s python style guide](https://google.github.io/styleguide/pyguide.html)
- python’s “official” style guide: [general style](https://peps.python.org/pep-0008/) and [docstrings](https://peps.python.org/pep-0257/)
- an example of enforcing style in a large, multi-language astro code project: [athena++ code](https://github.com/PrincetonUniversity/athena/wiki/Style-Guide)

# Basic styling

## Naming format

- variable and function: `lower_case_with_underscore`, but follow the convention in capitalizing variables; for example, `M_star`, `Sigma`
- class: `CamelCase` (and keep all cap acronyms capitalized, for example, `ALMAData`)

## Choosing a good name for your variable/function/class

- Use single letters (including Greek letters) and abbreviated subscripts if they are sufficiently commonly used (e.g., when most papers in your field use them).
Example: `a, e, Sigma, Omega, _init, _d, _in, _out`
- Otherwise, make sure the name is descriptive enough - don’t be afraid to use a long name.
Example: `minimum_disk_mass, outermost_semi_major_axis`
- If the name turns out to be extremely long (say, >4 words), try:
    - shortening words using common abbreviations: `minimum_disk_mass → M_d_min`
    - using a shorter name and add a comment: `initial_number_of_embryos_per_log_a → n_embryo_init # number of embryos per log a` (here `n` is commonly used for number density)
- Add a comment if there is any room of ambiguity; this could also be useful for labeling what unit system you use for this variable.
Example: `M_d_min = 0.01 # minimum disk mass in M_sun`
- Don’t shorten a name by omitting words.
Example: `outermost_semi_major_axis → outermost_axis`
- Don’t omit underscores.
Example: `a_out → aout`
- Don’t shorten Greek letters.
Example: `Sigma → Sig or S`
- Exception: Using a very short name is OK if it’s used repeatedly and only in a short block of code (e.g., when it will appear 20 times in the next 5 lines). But do add a comment to explain what this means when it first appears.
Example: `kp # kp = kappa = opacity`

## Comments

Commenting as you code is easier than adding comments at a later time. Add comments when

- You want to break a long script into steps:

```python
# step 1. do this
......
# step 2. do that
......
```

- The code is not self-explanatory, such as
    - a potentially ambiguous variable name
    - a reference for a equation you took from a paper
    - a very complex algorithm that does not seem intuitive even to yourself

## Docstrings for functions and classes

A docstring is a block of comment at the beginning of a function/class. In python, `"""..."""` is usually used for docstrings. As a rule of thumb, a docstring should allow you to understand what a function does and how to call it without reading the body of the function.

Always write the docstring before you write the body of the function. This forces yourself to have a clear idea of what you want to achieve before digging into the details. (Also, with a good docstring, AI can sometimes just write the function for you.)

What to include in a docstring:

- what the function/class does (can be omitted if the function/class name already says that very clearly)
- for function: inputs and outputs (including data types / units / dimensions of arrays)
- for class: attributes and (optionally) methods
- (optional) errors to be raised
- (optional) example use cases

Two examples from [my own code](https://github.com/wxu26/GIdisk2obs/blob/main/disk_model.py) (it doesn’t strictly follow these rules; but this gives you an idea of how much leeway one has in practice):

```python
def generate_disk_property_table(
    opacity_table, N_T_eff=30, N_Sigma_cs=30, T_eff_min=0.5,
    ):
    """
    Generate a table that maps T_eff and Sigma/cs to other local disk
    properties.

    Args:
      opacity_table: opacity table from generate_opacity_table()
      N_T_eff, N_Sigma_cs: resolution of the grid
      T_eff_min: minimum T_eff (max T_eff is determined automatically)

    Returns:
      a dictionary containing:
        T_eff_grid
        Sigma_cs_l, Sigma_cs_r: min/max Sigma/cs allowed at given T_eff
        x: normalized log(Sigma/cs) grid.
           we map x=[0,1] to the full range of allowed log(Sigma/cs) 
           uniformly.
        Sigma, cs, tau_p_mid, tau_r_mid, T_mid: local disk properties
           for given (T,x)
    """
```

```python
class DiskModel:
    """
    Parametrized disk model for generating radial porfiles of disk
    properties and flux density at given wavelengths.

    Attributes:
      (all in cgs)
      M: total mass
      Mstar: stellar mass
      Mdot: accretion rate
      Rd: disk size
      Q: Toomre Q
      
      R: radius grid (R[0]=0)
      Sigma, T_mid, tau_p_mid, tau_r_mid: radial profile at R[1:]
      MR: M(<R) profile at R[1:]
    """
```

For more examples, see section 3.8 of [google’s style guide](https://google.github.io/styleguide/pyguide.html).

## Function or class?

- If you see your functions frequently taking the same set of inputs, consider switching to a class where the functions are methods, or just store this set of inputs in a class / dictionary.
- If methods in a class barely needs to access attributes and other methods in the class, consider making it a function instead.

## Coding in jupyter notebook

During development, a notebook can look quite messy - that’s completely normal. But once in a while (say, when you know that you will need to pause the project for more than 24 hours) you should clean up the notebook to follow these standards:

- If you just rerun the whole notebook, it should give you the desired output. Don’t have earlier cells depending on results from later cells.
- Use markdown cells to (1) include a short description of what the whole notebook does and (2) divide the notebook into sections. Optionally, add longer descriptions and discussions, including instructions on how to modify the notebook to do slightly different things.
- Cell size: There is no universal guideline, but my recommendation is that each cell should contain a single class, a single function, or a single step (e.g., plotting a figure, or one step in a data processing pipeline). If what’s done in the cell is not obvious, add a comment (or a docstring, or a markdown cell) at the top of the cell explaining what the cell does.
- Notebook size: Generally speaking, a notebook does one “thing”: processing data, making a figure, making a set of figures, performing one certain analysis, etc. But what one considers as one “thing” can be somewhat flexible. If you find it difficult scrolling through it or running the whole notebook, it’s usually a sign that you want to break a notebook into multiple shorter ones.
- Naming notebooks: Just like functions, the names of notebooks need to mean something. (Don’t leave anything named untitled overnight…) There is no universal guideline for format (say, space vs underscore), so it’s fine as long as all your notebooks follow the same format.
- If a function will be used by more than one notebook, consider saving it in a script/module and import it.

# Files and folders

## Naming simulations and data files

We often have many simulations / data sets where each file/folder covers a different parameters. When naming these file/folders, it’s helpful to encode the parameters in the name. But unlike variable names where we want it to be as clear as possible, here we often want to prioritize conciseness since we don’t want the name get too long.

Example: `sim_a1_t2_w3` where a stands for initial a in au, t stands for duration in Myr, and w stands for width of the system (a_out/a_in)

When doing this kind of naming, include a readme file and document what the abbreviations mean.

## Folder structure

Some rules of thumb:

- If there are more than 10~20 different files in a folder, consider splitting it into different subfolders.
    - The number can be relaxed into ~1000 for folders containing data of the same type/format; e.g., a large batch of simulations.
- Splitting folders and organizing a project:
    - put data files and analysis scripts in different folders, when there are multiple of both
    - when there are more than a few .py module files, consider putting python modules and jupyter botebooks in different folders
- Deleting obsolete files:
    - As long as there is space available, don’t delete anything unless you are sure that they are extremely wrong and can’t possibly be useful anymore. Instead of deleting, just stash everything obsolete into a folder (and maybe add a readme on why they are obsolete).
    - A good time to delete obsolete files is when a project has been finalized and published.
- Temporary files and exploratory analysis: It’s OK to have temporary files, but do label them clearly and clean them regularly. For example, I often keep a scratch.ipynb for exploratory analysis, but the content is regularly deleted or moved to a more permanent notebook.

If you find it difficult to find an optimal way to structure your files and folders, that’s completely normal - this is a hard task. You can always ask your friends/mentors for suggestions.

# Sharing code and collaborating

## Minimal reproducible examples (MRE)

If you see an error in your code and want to ask others how to fix it, it’s useful to try making a minimal reproducible example where the error can be reproduces from a small, isolated script. This will help others (and yourself) isolate the cause of the problem and speed up debugging.

See a longer discussion on how to produce a MRE in [this article](https://stackoverflow.com/help/minimal-reproducible-example).

Sometimes it’s hard to produce a MRE that’s small enough to be shareable; maybe the code just has too many dependencies. But at least make sure that your error is reproducible.

## Readme

Think about it as a project-level docstring. Include things like

- What this project does
- How are folders and files organized
- How to reproduce this project, if it requires anything more than hitting run in your jupyter notebooks. (e.g., anything not in the repository that one needs to download? any special package dependencies?)
- References (related publications, who contributed, etc.)

## Version control

- use git
- commit often
- write informative commit messages

## Archiving projects

Github is usually enough for archiving analysis scripts and plots. For slightly larger datasets (say, GB-TB level), consider [zenodo](https://zenodo.org/) or [dataverse](https://dataverse.harvard.edu/). For even larger ones, it’s often hard to find a standardized solution.

Theoretically one shouldn’t just let data lie in some university cluster and die. In reality people (including myself) do that a lot…
