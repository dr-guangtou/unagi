#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = '1905'

from . import hsc
from . import sky
from . import mask
from . import task
from . import query
from . import config
from . import plotting

__all__ = ["query", "hsc", "task", "config", "sky", "mask", "plotting"]


def unagi():
    """Show random video about unagi on Youtube."""
    video_list = [
        "1sqLCUuMMfo", "p4KFCAX6X4o", "N26pjkM_z4A", "dn88LiPOKMc",
        "XUzsBV1xPNI", "TVEU-Pfj5eA", "7hnYiT23AEU", "4Qav8bdjCeg",
        "iaTr4V17Dwg",
    ]

    youtube_url = "https://www.youtube.com/watch?v="
    youtube_suffix = "?rel=0&amp;controls=0&amp;showinfo=0"
    video_urls = [youtube_url + v + youtube_suffix for v in video_list]

    import numpy as np

    try:
        _ = get_ipython
        from IPython.display import YouTubeVideo
        return YouTubeVideo(np.random.choice(video_list), width=560, height=315)
    except NameError:
        import webbrowser
        url = np.random.choice(video_urls)
        webbrowser.open(url, new=0, autoraise=True)
