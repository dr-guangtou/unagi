#!/usr/bin/env python
# -*- coding: utf-8 -*-

from . import hsc
from . import sky
from . import mask
from . import task
from . import query
from . import config

__all__ = ["query", "hsc", "task", "config", "sky", "mask"]

__version__ = "0.1.3"
__name__ = 'unagi'


def unagi():
    """Show random video about unagi on Youtube."""
    video_list = [
        "1sqLCUuMMfo", "p4KFCAX6X4o", "N26pjkM_z4A", "dn88LiPOKMc",
        "XUzsBV1xPNI", "TVEU-Pfj5eA", "7hnYiT23AEU", "4Qav8bdjCeg",
        "iaTr4V17Dwg", "pcVcJZUX0Dc", "lckYtOBTt48", "dgYeXXGYWtg",
        "E_TuYhsjJ_Y", "XTWchA5WlGQ", "xwNUtMIvYO4", "tSloyrM_XOg",
        "-0NsqEbE57Q"
    ]

    youtube_url = "https://www.youtube.com/watch?v="
    youtube_suffix = "?rel=0&amp;controls=0&amp;showinfo=0"
    video_urls = [youtube_url + v + youtube_suffix for v in video_list]

    import numpy as np

    try:
        ip = get_ipython()
        if ip.has_trait('kernel'):
            from IPython.display import YouTubeVideo
            return YouTubeVideo(np.random.choice(video_list), width=560, height=315)
        else:
            import webbrowser
            url = np.random.choice(video_urls)
            webbrowser.open(url, new=0, autoraise=True)
    except NameError:
        import webbrowser
        url = np.random.choice(video_urls)
        webbrowser.open(url, new=0, autoraise=True)
