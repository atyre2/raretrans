library(dplyr)
library(googledrive)
files <- list.files(R.home("doc"), pattern="NEWS", full.names=TRUE)[1:3]
# doesn't work?
drive_upload(media=files)
drive_upload(media=as_dribble(files))
# workaround:
newsDrive <- list()
for (i in seq_along(files)){
  newsDrive[[i]] <- drive_upload(media=files[i])
}
news <- bind_rows(newsDrive) # great! this works why did mine break?

# doesn't work with multiple files?
# docs say "...identifies the file(s) ..."
drive_share(news, role="reader", type="anyone")
# this works
drive_share(news[1,], role="reader", type="anyone")
# and this works but only the first one shareable (not surprising)
drive_share_link(news)
# cleanup
drive_rm(news)
# that worked too!
