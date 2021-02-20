#!/bin/sh

GREP_OPTIONS=''

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
finish () {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username: " username
    username=${username:-wt0881}
    read -s -p "Password: " password
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n00e115.zip"
    echo
    exit 1
}

prompt_credentials
  detect_app_approval() {
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 2 --netrc-file "$netrc" https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n00e115.zip -w %{http_code} | tail  -1`
    if [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    status=$(curl -s -z "$(date)" -w %{http_code} https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n00e115.zip | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

fetch_urls() {
  if command -v curl >/dev/null 2>&1; then
      setup_auth_curl
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -g -o $stripped_query_params -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  elif command -v wget >/dev/null 2>&1; then
      # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
      echo
      echo "WARNING: Can't find curl, use wget instead."
      echo "WARNING: Script may not correctly identify Earthdata Login integrations."
      echo
      setup_auth_wget
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document $stripped_query_params --keep-session-cookies -- $line && echo || exit_with_error "Command failed with error. Please retrieve the data manually."
      done;
  else
      exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
  fi
}

fetch_urls <<'EDSCEOF'
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w127.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n40w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w126.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w127.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w126.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n41w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w128.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w115.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w127.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w128.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w126.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w122.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w128.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w119.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n45w120.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n44w125.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n49w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w123.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n47w118.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n48w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n42w121.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n43w116.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n46w117.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n50w126.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w124.zip
https://e4ftl01.cr.usgs.gov//DP132/MEASURES/NASADEM_HGT.001/2000.02.11/NASADEM_HGT_n51w119.zip
EDSCEOF
