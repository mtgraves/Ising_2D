import pylab as pl

# =============================================================================
def darkPlots(s):
    """
    Sets all parts of a standard x vs. y plot
    to white, including spines, tick marks, labels,
    title, and axis values
    """
    return (s.patch.set_facecolor('none'),
            s.spines['bottom'].set_color('white'), 
            s.spines['top'].set_color('white'),
            s.spines['left'].set_color('white'),
            s.spines['right'].set_color('white'),
            s.tick_params(axis='x', color='white', labelcolor='white'),
            s.tick_params(axis='y', color='white', labelcolor='white'),
            s.title.set_color('white'),
            s.xaxis.label.set_color('white'),
            s.yaxis.label.set_color('white'))
